
import torch
from torch import nn
import torchvision

NUM_INPUT_CHANNELS = 3

# Define model
class HrdModel(nn.Module):
    def __init__(self, dropout_rate, num_tumor_types):
        super().__init__()

        self.resnet = torchvision.models.resnet18(weights=torchvision.models.ResNet18_Weights.DEFAULT)

        with torch.no_grad():
            # add dropout before each relu
            print(f"adding dropout after each relu layer, rate={dropout_rate}")
            append_dropout(self.resnet, dropout_rate)

            # change to 4 channels
            # https://discuss.pytorch.org/t/how-to-modify-the-input-channels-of-a-resnet-model/2623/8

            new_in_channels = NUM_INPUT_CHANNELS
            layer = self.resnet.conv1

            # Creating new Conv2d layer
            new_layer = nn.Conv2d(in_channels=new_in_channels,
                                  out_channels=layer.out_channels,
                                  kernel_size=layer.kernel_size,
                                  stride=layer.stride,
                                  padding=layer.padding,
                                  bias=layer.bias)

            copy_weights = 0  # Here will initialize the weights from new channel with the red channel weights

            # Copying the weights from the old to the new layer
            new_layer.weight[:, :layer.in_channels, :, :] = layer.weight.clone()

            # Copying the weights of the `copy_weights` channel of the old layer to the extra channels of the new layer
            for i in range(new_in_channels - layer.in_channels):
                channel = layer.in_channels + i
                new_layer.weight[:, channel:channel + 1, :, :] = layer.weight[:, copy_weights:copy_weights + 1,
                                                                 ::].clone()
            new_layer.weight = nn.Parameter(new_layer.weight)

            self.resnet.conv1 = new_layer
            
            self.resnet.fc = nn.Identity(512)

            # now make fc layers
            self.fc_stack = nn.Sequential(
                nn.Linear(512 + num_tumor_types, 64),
                nn.ReLU(),
                nn.BatchNorm1d(64),
                nn.Dropout(dropout_rate),
                nn.Linear(64, 32),
                nn.ReLU(),
                nn.BatchNorm1d(32),
                nn.Dropout(dropout_rate),
                nn.Linear(32, 16),
                nn.ReLU(),
                nn.BatchNorm1d(16),
                nn.Dropout(dropout_rate),
                nn.Linear(16, 1))

    def forward(self, image_x, linear_x):
        
        #print(seq_x.shape)
        x = self.resnet(image_x)

        #print(x.shape)
        
        # concat the sequential input with the linear data
        x = torch.cat((x, linear_x), dim=1)
        
        # uncomment follow to work out how big the first linear layer needs to be
        # print(x.shape)
        
        x = self.fc_stack(x)
        return x

def append_dropout(model, rate):
    for name, module in model.named_children():
        if len(list(module.children())) > 0:
            append_dropout(module, rate)
        if isinstance(module, nn.ReLU):
            new = nn.Sequential(module, nn.Dropout2d(p=rate))
            setattr(model, name, new)
