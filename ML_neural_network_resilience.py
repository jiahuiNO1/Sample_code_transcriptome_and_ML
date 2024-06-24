#Binary task adversarial neural network for resilience

import torch
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from torch.utils.data import DataLoader, TensorDataset
import torch.nn as nn
import torch.optim as optim
from sklearn.metrics import roc_auc_score
import optuna
from torch.autograd import Function
from captum.attr import IntegratedGradients

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
# device = "cpu"
print("Using device:", device)

# Define the paths to your train files
train_path_x1 = "D:/AD_resilience_deep_learning/resilience_ML_codes/simulated_data/matched_simulated_low_risk_set.train.txt"
train_path_x2 = "D:/AD_resilience_deep_learning/resilience_ML_codes/simulated_data/matched_simulated_resilient_set.train.txt"
valid_path_x1 = "D:/AD_resilience_deep_learning/resilience_ML_codes/simulated_data/matched_simulated_low_risk_set.valid.txt"
valid_path_x2 = "D:/AD_resilience_deep_learning/resilience_ML_codes/simulated_data/matched_simulated_resilient_set.valid.txt"

# Load the training datasets
train_ds_x1 = pd.read_csv(train_path_x1)
train_ds_x2 = pd.read_csv(train_path_x2)

# Concatenate along the sample axis assuming they have the same features
train_ds = pd.concat([train_ds_x1, train_ds_x2], axis=0, ignore_index=True) #rbind

# Extract features (x1 and x2) and labels (y)
x1_train, x2_train, y_train_x1, y_train_x2 = (
    train_ds_x1.values[:, 5:-2],
    train_ds_x2.values[:, 5:-2],
    train_ds_x1.values[:, 1],
    train_ds_x2.values[:, 1]
)

# Transform into NumPy array
x1_train, x2_train = np.asarray(x1_train).astype('float32'), np.asarray(x2_train).astype('float32')
# Encode categorical (textual) data into numerical labels
y_train_x1, y_train_x2 = LabelEncoder().fit_transform(y_train_x1), LabelEncoder().fit_transform(y_train_x2)

# Load the validation datasets (similarly as for training)
valid_ds_x1 = pd.read_csv(valid_path_x1)
valid_ds_x2 = pd.read_csv(valid_path_x2)

# Concatenate along the sample axis for validation data
valid_ds = pd.concat([valid_ds_x1, valid_ds_x2], axis=0, ignore_index=True)

# Extract features (x1 and x2) and labels (y)
x1_valid, x2_valid, y_valid_x1, y_valid_x2 = (
    valid_ds_x1.values[:, 5:-2],
    valid_ds_x2.values[:, 5:-2],
    valid_ds_x1.values[:, 1],
    valid_ds_x2.values[:, 1]
)

x1_valid, x2_valid = np.asarray(x1_valid).astype('float32'), np.asarray(x2_valid).astype('float32')
y_valid_x1, y_valid_x2 = LabelEncoder().fit_transform(y_valid_x1), LabelEncoder().fit_transform(y_valid_x2)

# Print the shapes of the data
print("x1 train shape:", x1_train.shape)
print("x2 train shape:", x2_train.shape)
print("y1 train shape:", y_train_x1.shape)
print("y2 train shape:", y_train_x2.shape)

# Create DataLoader objects for each subset
batch_size = 2056  # Adjust the batch size as needed

# Training DataLoader for x1
# -creates a PyTorch Dataset object (TensorDataset) containing a tensor representing the input data (x1_train).
# -a tensor is a multi-dimensional array or matrix
# -a 3-dimensional tensor with shape (2, 3, 4) has 2 layers, each containing a 3x4 matrix.
train_dataset_x1 = TensorDataset(torch.tensor(x1_train, dtype=torch.float32).to(device), torch.tensor(y_train_x1, dtype=torch.float32).to(device))
# Load and iterate over data in batches, randomly shuffles the data at the start of each epoch 
train_loader_x1 = DataLoader(train_dataset_x1, batch_size=batch_size, shuffle=True)

# Training DataLoader for x2
train_dataset_x2 = TensorDataset(torch.tensor(x2_train, dtype=torch.float32).to(device), torch.tensor(y_train_x2, dtype=torch.float32).to(device))
train_loader_x2 = DataLoader(train_dataset_x2, batch_size=batch_size, shuffle=True)

# Validation DataLoader for x1
val_dataset_x1 = TensorDataset(torch.tensor(x1_valid, dtype=torch.float32).to(device), torch.tensor(y_valid_x1, dtype=torch.float32).to(device))
val_loader_x1 = DataLoader(val_dataset_x1, batch_size=batch_size, shuffle=False)

# Validation DataLoader for x2
val_dataset_x2 = TensorDataset(torch.tensor(x2_valid, dtype=torch.float32).to(device), torch.tensor(y_valid_x2, dtype=torch.float32).to(device))
val_loader_x2 = DataLoader(val_dataset_x2, batch_size=batch_size, shuffle=False)

#
class ReverseLayerF(Function):
    @staticmethod
    def forward(ctx, x, gamma):
        ctx.gamma = gamma  # Store gamma in the context
        return x.view_as(x)

    @staticmethod
    def backward(ctx, grad_output):
        output = grad_output.neg() * ctx.gamma  # Reverse the gradients with respect to gamma
        return output, None


# Define binary classification model
class BinaryClassifier(nn.Module):
    def __init__(self, input_dim_x1, input_dim_x2, num_layers, layer_size, dropout_rate, gamma):
        super(BinaryClassifier, self).__init__()

        # Shared layers
        shared_layers = []
        shared_layers.append(nn.Linear(input_dim_x1, layer_size))
        shared_layers.append(nn.ReLU())

        for _ in range(num_layers - 1):
            shared_layers.append(nn.Linear(layer_size, layer_size))
            shared_layers.append(nn.ReLU())
            shared_layers.append(nn.Dropout(p=dropout_rate))  # Add a dropout layer after each ReLU

        self.shared_layers = nn.Sequential(*shared_layers)

        # Task-specific layers for x1
        self.task_x1 = nn.Linear(layer_size, 1)
        self.sigmoid_x1 = nn.Sigmoid()

        # Task-specific layers for x2
        self.task_x2 = nn.Linear(layer_size, 1)
        self.sigmoid_x2 = nn.Sigmoid()

        # Add the ReverseGrad layer with the specified gamma
        self.reverse_grad = ReverseLayerF.apply
        self.gamma = gamma

        self.to(device)  # Move the entire model to GPU if available

    def forward(self, x1, x2):
        # Shared layers
        x1_shared = self.shared_layers(x1)
        x2_shared = self.shared_layers(x2)

        # Task-specific layers for x1
        # includes reverse layer
        x1_rev = self.reverse_grad(x1_shared, self.gamma)
        x1_output = self.sigmoid_x1(self.task_x1(x1_rev))

        # Task-specific layers for x2
        x2_output = self.sigmoid_x2(self.task_x2(x2_shared))

        return x1_output.to(device), x2_output.to(device)  # Move the outputs to GPU

from sklearn.metrics import accuracy_score

#
def evaluate_binary_classifier(model, dataloader_x1, dataloader_x2, criterion):
    model.eval()
    total_loss = 0.0
    #correct_x1 = 0
    #correct_x2 = 0
    #total = 0
    predicted_probabilities_x1 = []  # To store predicted probabilities for x1
    predicted_probabilities_x2 = []  # To store predicted probabilities for x2
    true_labels_x1 = []  # To store true labels for x1
    true_labels_x2 = []  # To store true labels for x2

    with torch.no_grad():
        for (inputs_x1, labels_x1), (inputs_x2, labels_x2) in zip(dataloader_x1, dataloader_x2):
            
            outputs_x1, outputs_x2 = model(inputs_x1, inputs_x2)

            loss_x1 = criterion(outputs_x1.squeeze(), labels_x1)
            loss_x2 = criterion(outputs_x2.squeeze(), labels_x2)

            # Convert predictions to binary (0 or 1) based on a threshold (e.g., 0.5)
            #predictions_x1 = (outputs_x1 >= 0.5).float()
            #predictions_x2 = (outputs_x2 >= 0.5).float()

            total_loss += (loss_x1.item() + loss_x2.item()) / 2.0

            predicted_probabilities_x1.extend(outputs_x1.tolist())  # Append predicted probabilities for x1
            predicted_probabilities_x2.extend(outputs_x2.tolist())  # Append predicted probabilities for x2

            true_labels_x1.extend(labels_x1.tolist())  # Append true labels for x1
            true_labels_x2.extend(labels_x2.tolist())  # Append true labels for x2

    average_loss = total_loss / len(dataloader_x1)  # Use the length of dataloader_x1 for consistency

    # Calculate AUC
    auc_x1 = roc_auc_score(true_labels_x1, predicted_probabilities_x1)
    auc_x2 = roc_auc_score(true_labels_x2, predicted_probabilities_x2)

    # Calculate sklearn accuracy using all predictions
    accuracy_x1 = accuracy_score(true_labels_x1, (np.array(predicted_probabilities_x1) >= 0.5).astype(float))
    accuracy_x2 = accuracy_score(true_labels_x2, (np.array(predicted_probabilities_x2) >= 0.5).astype(float))

    return accuracy_x1, accuracy_x2, average_loss, auc_x1, auc_x2


# Function to evaluate the single_task_model
def evaluate_single_task_model(model, dataloader_x2, criterion):
    model.eval()
    total_loss = 0.0
    #correct_x2 = 0
    #total = 0
    predicted_probabilities_x2 = []  # To store predicted probabilities for x2
    true_labels_x2 = []  # To store true labels for x2

    with torch.no_grad():
        for inputs_x2, labels_x2 in dataloader_x2:
            outputs_x2 = model(inputs_x2)

            loss_x2 = criterion(outputs_x2.squeeze(), labels_x2)

            # Convert predictions to binary (0 or 1) based on a threshold (e.g., 0.5)
            #predictions_x2 = (outputs_x2 >= 0.5).float()

            total_loss += loss_x2.item()

            predicted_probabilities_x2.extend(outputs_x2.tolist())  # Append predicted probabilities for x2
            true_labels_x2.extend(labels_x2.tolist())  # Append true labels for x2

    accuracy_x2 = accuracy_score(true_labels_x2, (np.array(predicted_probabilities_x2) >= 0.5).astype(float))

    average_loss = total_loss / len(dataloader_x2)

    auc_x2 = roc_auc_score(true_labels_x2, predicted_probabilities_x2)

    return accuracy_x2, average_loss, auc_x2

#
def objective(trial):
    config = {
        "lr": trial.suggest_float("lr", 1e-5, 1e-1, log=True),
        "epochs": trial.suggest_int("epochs", 10, 1000),
        "num_layers": trial.suggest_int("num_layers", 1, 5),
        "layer_size": trial.suggest_int("layer_size", 32, 5000),
        "dropout_rate": trial.suggest_float("dropout_rate", 0.0, 0.7),
        # set reverse layer parameter
        "gamma": trial.suggest_float("gamma", 0.0, 3.0),  # Optimize gamma as a parameter
    }

    print(f"Hyperparameters: {config}")

    # Initialize a list to store AUC values across multiple model initializations
    auc_values = []
    auc_values2 = []
    train_auc_values = []
    train_auc_values2 = []

    # Run multiple trials to ensure consistency across model initializations
    for trial_num in range(3):  # You can adjust the number of trials

        # Move entire model and associated tensors to GPU
        model = BinaryClassifier(x1_train.shape[1], x2_train.shape[1], config["num_layers"], config["layer_size"],
                                 config["dropout_rate"], config["gamma"])
        model = model.to(device)

        optimizer = optim.SGD(model.parameters(), lr=config["lr"])
        criterion = nn.BCELoss() #Binary Cross Entropy Loss

        for epoch in range(config["epochs"]):
            model.train()

            # Iterate over both loaders simultaneously
            for (inputs_x1, labels_x1), (inputs_x2, labels_x2) in zip(train_loader_x1, train_loader_x2):
                # Ensure labels are 1D
                labels_x1 = labels_x1.view(-1)
                labels_x2 = labels_x2.view(-1)

                optimizer.zero_grad()
                outputs_x1, outputs_x2 = model(inputs_x1, inputs_x2)
                loss_x1 = criterion(outputs_x1.squeeze(), labels_x1)
                loss_x2 = criterion(outputs_x2.squeeze(), labels_x2)
                loss = (loss_x1 + loss_x2) / 2.0
                loss.backward() #compute gradient of the loss
                optimizer.step() #update weights 

            # Evaluation on the training set
            model.eval()
            train_accuracy_x1, train_accuracy_x2, train_loss, train_auc_x1, train_auc_x2 = evaluate_binary_classifier(model, train_loader_x1, train_loader_x2, criterion)

            # Evaluation on the validation set
            model.eval()
            val_accuracy_x1, val_accuracy_x2, val_loss, val_auc_x1, val_auc_x2 = evaluate_binary_classifier(model, val_loader_x1, val_loader_x2, criterion)

            # Report the AUC value for x2 for each trial
            train_auc_values.append(train_auc_x2)
            train_auc_values2.append(train_auc_x1)
            auc_values.append(val_auc_x2)
            auc_values2.append(val_auc_x1)



    # Calculate the mean AUC value across multiple trials
    train_mean_auc = sum(train_auc_values) / len(train_auc_values)
    print(f"Mean train AUC main task: {train_mean_auc}")
    train_mean_auc2 = sum(train_auc_values2) / len(train_auc_values2)
    print(f"Mean train AUC adversarial task: {train_mean_auc2}")
    mean_auc = sum(auc_values) / len(auc_values)
    print(f"Mean valid AUC main task: {mean_auc}")
    mean_auc2 = sum(auc_values2) / len(auc_values2)
    print(f"Mean valid AUC adversarial task: {mean_auc2}")
    # Calculate difference between mean AUCs
    mean_auc_difference = mean_auc - (abs(0.5 - mean_auc2)+0.5)

    print(f"Mean valid AUC difference across trials: {mean_auc_difference}")
    return mean_auc_difference


##
study_name = "mt_opt_resilience_adv_tuning_4.11.24"
custom_storage_location = "sqlite:///D:\\AD_resilience_deep_learning\\resilience_ML_codes\\results\\optuna.db"


study = optuna.create_study(
    direction='maximize',
    study_name=study_name,
    storage=custom_storage_location,
    load_if_exists=True
)


# Continue the optimization with additional trials
study.optimize(objective, n_trials=0)  # You can adjust the number of additional trials

# Access the best parameters
best_params = study.best_params
print("Best Parameters:", best_params)


# Train the best model on the entire training dataset
# You might need to adjust the number of epochs and other training parameters
best_model = BinaryClassifier(x1_train.shape[1], x2_train.shape[1], best_params["num_layers"], best_params["layer_size"], best_params["dropout_rate"], best_params["gamma"])
best_model = best_model.to(device)

best_optimizer = optim.SGD(best_model.parameters(), lr=best_params["lr"])
best_criterion = nn.BCELoss() #Binary Cross Entropy Loss

for epoch in range(best_params["epochs"]):
    best_model.train()

    # Iterate over both loaders simultaneously
    for (inputs_x1, labels_x1), (inputs_x2, labels_x2) in zip(train_loader_x1, train_loader_x2):
        # Ensure labels are 1D
        labels_x1 = labels_x1.view(-1)
        labels_x2 = labels_x2.view(-1)

        best_optimizer.zero_grad()
        outputs_x1, outputs_x2 = best_model(inputs_x1, inputs_x2)
        loss_x1 = best_criterion(outputs_x1.squeeze(), labels_x1)
        loss_x2 = best_criterion(outputs_x2.squeeze(), labels_x2)
        loss = (loss_x1 + loss_x2) / 2.0
        loss.backward()
        best_optimizer.step()

# Evaluate the best model on the validation set
best_model.eval()
best_val_accuracy_x1, best_val_accuracy_x2,best_val_loss, best_val_auc_x1, best_val_auc_x2 = evaluate_binary_classifier(best_model, val_loader_x1, val_loader_x2, best_criterion)

# Report the performance of the best model
print("Best Model Performance:")
print("Validation Accuracy (x1): {:.2f}%".format(best_val_accuracy_x1))
print("Validation Accuracy (x2): {:.2f}%".format(best_val_accuracy_x2))
print("Validation Loss: {:.4f}".format(best_val_loss))
print("Validation AUC (x1): {:.4f}".format(best_val_auc_x1))
print("Validation AUC (x2): {:.4f}".format(best_val_auc_x2))

# Save the weights of the best_model
torch.save(best_model.state_dict(), 'best_model_weights_adv.pth') #retrieve the current state of the model parameters


# Define a new model that only uses x2 input and task
class SingleTaskModel(nn.Module):
    def __init__(self, input_dim, num_layers, layer_size, dropout_rate):
        super(SingleTaskModel, self).__init__()

        # Shared layers
        shared_layers = []
        shared_layers.append(nn.Linear(input_dim, layer_size))
        shared_layers.append(nn.ReLU())

        for _ in range(num_layers - 1):
            shared_layers.append(nn.Linear(layer_size, layer_size))
            shared_layers.append(nn.ReLU())
            shared_layers.append(nn.Dropout(p=dropout_rate))  # Add a dropout layer after each ReLU

        self.shared_layers = nn.Sequential(*shared_layers)

        # Task-specific layer for x2
        self.task_x2 = nn.Linear(layer_size, 1)
        self.sigmoid_x2 = nn.Sigmoid()

        self.to(device)

    def forward(self, x2):
        # Shared layers
        x2_shared = self.shared_layers(x2)

        # Task-specific layer for x2
        x2_output = self.sigmoid_x2(self.task_x2(x2_shared))

        return x2_output.to(device)

# Instantiate the new model
single_task_model = SingleTaskModel(x2_train.shape[1], best_params["num_layers"], best_params["layer_size"], best_params["dropout_rate"])

# Load the saved weights into the new model, skipping keys not present in the new model
state_dict = torch.load('best_model_weights_adv.pth')
model_dict = single_task_model.state_dict()

# Filter out unnecessary keys
state_dict = {k: v for k, v in state_dict.items() if k in model_dict}

# Update the model's state_dict
single_task_model.load_state_dict(state_dict)

# Evaluate the single_task_model on the validation set
single_task_model.eval()
single_task_val_accuracy_x2, single_task_val_loss, single_task_val_auc_x2 = evaluate_single_task_model(single_task_model, val_loader_x2, best_criterion)

# Report the performance of the single_task_model
print("Single Task Model Performance:")
print("Validation Accuracy (x2): {:.2f}%".format(single_task_val_accuracy_x2))
print("Validation Loss: {:.4f}".format(single_task_val_loss))
print("Validation AUC (x2): {:.4f}".format(single_task_val_auc_x2))


#
def get_integrated_gradients_single_task(model, inputs_x2, baseline_x2, n_steps=50, device=device):
    model.eval()

    # Ensure the model and inputs are on the same device
    model = model.to(device)
    inputs_x2 = inputs_x2.to(device)
    baseline_x2 = baseline_x2.to(device)

    # Initialize IntegratedGradients
    ig = IntegratedGradients(model) #quantifies the importance of each feature by integrating the gradients of the model's prediction with respect to the input features over a path from a baseline to the input.

    # Compute attributions for x2
    attributions_x2 = ig.attribute(inputs_x2, baselines=baseline_x2, n_steps=n_steps)

    return attributions_x2


# Assuming you have a DataFrame `valid_ds_x2` with feature names
feature_names = valid_ds_x2.columns[5:-2]  # Adjust the column index based on your data

# Loop through the entire validation dataset
attributions_cases_list = []
attributions_controls_list = []

for i in range(len(valid_ds_x2)):
    # Choose a sample from the validation set for analysis
    x2_sample = x2_valid[i:i + 1]
    baseline_x2 = torch.zeros_like(torch.tensor(x2_sample)).to(device)
    x2_tensor = torch.tensor(x2_sample, dtype=torch.float32).to(device)

    # Perform feature importance analysis
    attributions_x2 = get_integrated_gradients_single_task(single_task_model, x2_tensor, baseline_x2, n_steps=50)

    # Get the label for the current sample
    label = valid_ds_x2.loc[i, 'PHENO']  # Replace 'Your_Label_Column_Name' with the actual column name

    # Move the tensor to the CPU before converting to NumPy array
    attributions_x2_np = attributions_x2.cpu().numpy()

    # Append attributions to the corresponding list
    if label == 1:  # Assuming 1 is the label for cases
        attributions_cases_list.append(attributions_x2_np)
    elif label == 0:  # Assuming 0 is the label for controls
        attributions_controls_list.append(attributions_x2_np)

# Concatenate the lists into DataFrames
df_attributions_cases = pd.DataFrame(np.concatenate(attributions_cases_list, axis=0), columns=feature_names)
df_attributions_controls = pd.DataFrame(np.concatenate(attributions_controls_list, axis=0), columns=feature_names)

# Save DataFrames to CSV files
df_attributions_cases.to_csv('attributions_cases_adv_4.11.csv', index=False)
df_attributions_controls.to_csv('attributions_controls_adv_4.11.csv', index=False)



### do it again for other input
# Loop through the entire validation dataset
attributions_cases_list = []
attributions_controls_list = []

for i in range(len(valid_ds_x1)):
    # Choose a sample from the validation set for analysis
    x1_sample = x1_valid[i:i + 1]
    baseline_x1 = torch.zeros_like(torch.tensor(x1_sample)).to(device)
    x1_tensor = torch.tensor(x1_sample, dtype=torch.float32).to(device)

    # Perform feature importance analysis
    attributions_x1 = get_integrated_gradients_single_task(single_task_model, x1_tensor, baseline_x1, n_steps=50)

    # Get the label for the current sample
    label = valid_ds_x1.loc[i, 'PHENO']  # Replace 'Your_Label_Column_Name' with the actual column name

    # Move the tensor to the CPU before converting to NumPy array
    attributions_x1_np = attributions_x1.cpu().numpy()

    # Append attributions to the corresponding list
    if label == 1:  # Assuming 1 is the label for cases
        attributions_cases_list.append(attributions_x1_np)
    elif label == 0:  # Assuming 0 is the label for controls
        attributions_controls_list.append(attributions_x1_np)

# Concatenate the lists into DataFrames
df_attributions_cases = pd.DataFrame(np.concatenate(attributions_cases_list, axis=0), columns=feature_names)
df_attributions_controls = pd.DataFrame(np.concatenate(attributions_controls_list, axis=0), columns=feature_names)

# Save DataFrames to CSV files
df_attributions_cases.to_csv('attributions_cases_adv_secondary_4.11.csv', index=False)
df_attributions_controls.to_csv('attributions_controls_adv_secondary_4.11.csv', index=False)


