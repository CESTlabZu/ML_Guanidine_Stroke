import torch
import torch.nn as nn

class StrokeNN(nn.Module):
    def __init__(self, input_dim):
        super(StrokeNN, self).__init__()
        # Define fully connected layers with batch normalization
        self.fc1 = nn.Linear(input_dim, 512)  # Input layer
        self.fc2 = nn.Linear(512, 64)
        self.fc3 = nn.Linear(64, 64)
        self.fc4 = nn.Linear(64, 16)
        self.output = nn.Linear(16, 2)  # Output layer

    def forward(self, x):
        x = torch.relu(self.fc1(x))
        x = torch.relu(self.fc2(x))
        x = torch.relu(self.fc3(x))
        x = torch.relu(self.fc4(x))
        x = self.output(x)

        return x
