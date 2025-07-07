import torch
from torch import nn
import torchvision
from enum import Enum

IMAGE_SIZE = 512
NUM_INPUT_CHANNELS = 3
NUM_DRIVERS = 729
NUM_DRIVERS = 671
NUM_SNVS = 8
#NUM_SNVS = 105

cancer_types = ['Colorectum/Small intestine/Appendix', 'Breast: Other', 'Prostate',
                'Skin: Melanoma', 'Urothelial tract', 'Lung: Non-small cell: LUAD',
                'Gynecologic: Ovary/Fallopian tube', 'Esophagus/Stomach',
                'HPB: Pancreas', 'Bone/Soft tissue: Other', 'Kidney: Other',
                'Breast: Triple negative', 'Anogenital',
                'HPB: Bile duct/Gallbladder', 'CNS: Glioma', 'Mesothelium',
                'Bone/Soft tissue: GIST', 'Bone/Soft tissue: Leiomyosarcoma',
                'NET: Colorectum/Small intestine', 'Lung: Small cell',
                'HPB: Liver', 'Head and neck: Other', 'Skin: Other', 'NET: Lung',
                'Gynecologic: Endometrium', 'Bone/Soft tissue: Liposarcoma',
                'Bone/Soft tissue: Undiff. sarcoma', 'NET: Pancreas',
                'Lung: Non-small cell: LUSC', 'Thyroid gland',
                'Head and neck: Salivary gland', 'Lymphoid tissue',
                'Bone/Soft tissue: Osteosarcoma', 'Head and neck: Adenoid cystic',
                'CNS: Medulloblastoma', 'Kidney: Chromophobe'][:20]


# Define model
class VCuppaModel(nn.Module):
    def __init__(self, cnn_dropout_rate, linear_dropout_rate):
        super().__init__()

        self.resnet = torchvision.models.resnet18(weights=torchvision.models.ResNet18_Weights.DEFAULT)

        with torch.no_grad():

            # add dropout before each relu
            print(f"adding dropout after each relu layer, rate={cnn_dropout_rate}")
            append_dropout(self.resnet, cnn_dropout_rate)

            # adding batchnorm after image helps boost performance
            self.resnet.conv1 = nn.Sequential(
                nn.BatchNorm2d(3),  # 3 channels
                self.resnet.conv1)

            self.resnet.fc = nn.Identity(512)

            # now make fc layers
            self.fc_stack = nn.Sequential(
                nn.Linear(512 + NUM_DRIVERS + NUM_SNVS, 256),
                nn.ReLU(),
                nn.BatchNorm1d(256),
                nn.Dropout(linear_dropout_rate),
                nn.Linear(256, 128),
                nn.ReLU(),
                nn.BatchNorm1d(128),
                nn.Dropout(linear_dropout_rate),
                nn.Linear(128, 64),
                nn.ReLU(),
                nn.BatchNorm1d(64),
                nn.Linear(64, len(cancer_types)))

    def forward(self, image_x, drivers_x):
        # print(seq_x.shape)
        x = self.resnet(image_x)

        # print(x.shape)

        # concat the sequential input with the linear data
        x = torch.cat((x, drivers_x), dim=1)

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
