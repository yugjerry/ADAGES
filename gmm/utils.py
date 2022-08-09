{
  "nbformat": 4,
  "nbformat_minor": 0,
  "metadata": {
    "colab": {
      "name": "utils.ipynb",
      "provenance": [],
      "authorship_tag": "ABX9TyO185zQm/kE/UGJi4I52fNQ",
      "include_colab_link": true
    },
    "kernelspec": {
      "name": "python3",
      "display_name": "Python 3"
    },
    "language_info": {
      "name": "python"
    }
  },
  "cells": [
    {
      "cell_type": "markdown",
      "metadata": {
        "id": "view-in-github",
        "colab_type": "text"
      },
      "source": [
        "<a href=\"https://colab.research.google.com/github/yugjerry/ADAGES/blob/master/gmm/utils.py\" target=\"_parent\"><img src=\"https://colab.research.google.com/assets/colab-badge.svg\" alt=\"Open In Colab\"/></a>"
      ]
    },
    {
      "cell_type": "code",
      "execution_count": null,
      "metadata": {
        "id": "lDi2-hIxrMY7"
      },
      "outputs": [],
      "source": [
        "import numpy as np\n",
        "import numpy.linalg as LA\n",
        "import torch\n",
        "import torch.nn as nn\n",
        "from torch.utils.data import Dataset, DataLoader\n",
        "import torch.nn.functional as F\n",
        "from scipy import stats\n",
        "import matplotlib.pyplot as plt"
      ]
    },
    {
      "cell_type": "code",
      "source": [
        "## architecture\n",
        "class Projector(nn.Module): ## a linear projector\n",
        "    def __init__(self, input_dim, output_dim):\n",
        "        super(Projector, self).__init__()\n",
        "        self.input_dim = input_dim \n",
        "        self.output_dim = output_dim\n",
        "        self.linear = nn.Linear(self.input_dim, self.output_dim, bias=False)\n",
        "\n",
        "    def forward(self, x):\n",
        "        z = self.linear(x)\n",
        "        return z\n",
        "    \n",
        "\n",
        "class Projector2(nn.Module): ## 2-layer MLP projector\n",
        "    def __init__(self, input_dim, output_dim, hidden_dim):\n",
        "        super(Projector2, self).__init__()\n",
        "        self.input_dim = input_dim \n",
        "        self.hidden_dim = hidden_dim\n",
        "        self.output_dim = output_dim\n",
        "        self.linear1 = nn.Linear(self.input_dim, self.hidden_dim, bias=False)\n",
        "        self.linear = nn.Linear(self.hidden_dim, self.output_dim, bias=False)\n",
        "\n",
        "    def forward(self, x):\n",
        "        x = self.linear1(x)\n",
        "        z = self.linear(x)\n",
        "        return z\n",
        "    \n",
        "class Projector3(nn.Module): ## multilayer MLP projector\n",
        "    def __init__(self, input_dim, output_dim, hidden_dim):\n",
        "        super(Projector3, self).__init__()\n",
        "        self.input_dim = input_dim\n",
        "        self.output_dim = output_dim\n",
        "        self.hidden_dim = hidden_dim\n",
        "        current_dim = input_dim\n",
        "        self.linear = nn.Linear(self.hidden_dim[-1], self.output_dim, bias=False)\n",
        "        self.layers = nn.ModuleList()\n",
        "        for hdim in hidden_dim:\n",
        "            self.layers.append(nn.Linear(current_dim, hdim, bias=False))\n",
        "            current_dim = hdim\n",
        "        self.layers.append(self.linear)\n",
        "    def forward(self, x):\n",
        "        for layer in self.layers[:-1]:\n",
        "            x = layer(x)\n",
        "        out = self.linear(x)\n",
        "        return out "
      ],
      "metadata": {
        "id": "hNZ1MzV6rTZX"
      },
      "execution_count": null,
      "outputs": []
    },
    {
      "cell_type": "code",
      "source": [
        "def info_nce_loss(features, temperature, batch_size): # copies the code from SimCLR\n",
        "    labels = torch.cat([torch.arange(batch_size) for i in range(num_aug)], dim=0)\n",
        "    labels = (labels.unsqueeze(0) == labels.unsqueeze(1)).float()\n",
        "    features = F.normalize(features, dim=1)\n",
        "    similarity_matrix = torch.matmul(features, features.T)\n",
        "    #assert similarity_matrix.shape == (num_aug * batch_size, num_aug * batch_size)\n",
        "    #assert similarity_matrix.shape == labels.shape\n",
        "    mask = torch.eye(labels.shape[0], dtype=torch.bool)\n",
        "    labels = labels[~mask].view(labels.shape[0], -1)\n",
        "    similarity_matrix = similarity_matrix[~mask].view(similarity_matrix.shape[0], -1)\n",
        "    #assert similarity_matrix.shape == labels.shape\n",
        "\n",
        "    # select and combine multiple positives\n",
        "    positives = similarity_matrix[labels.bool()].view(labels.shape[0], -1)\n",
        "\n",
        "    # select only the negatives the negatives\n",
        "    negatives = similarity_matrix[~labels.bool()].view(similarity_matrix.shape[0], -1)\n",
        "\n",
        "    logits = torch.cat([positives, negatives], dim=1)\n",
        "    labels = torch.zeros(logits.shape[0], dtype=torch.long)\n",
        "    logits = logits / temperature\n",
        "\n",
        "    return logits, labels"
      ],
      "metadata": {
        "id": "ea3hpYS0rYiw"
      },
      "execution_count": null,
      "outputs": []
    },
    {
      "cell_type": "code",
      "source": [
        "## evaluation\n",
        "def calc_linear_Guassian_prob(mu, Sigma, a, b, side='positive'):\n",
        "    # Calculates the probability of P(<a,x> + b >=0) where x follows N(mu, Sigma) if side == 'positive\n",
        "    # P(<a,x> + b <= 0) if side == 'negative'\n",
        "    assert (mu.shape[0] == Sigma.shape[0]) and (mu.shape[0] == a.shape[0])\n",
        "    U, s, Vh = np.linalg.svd(Sigma)\n",
        "    assert np.all(s > 0)\n",
        "    Sigma_sqrt = np.matmul(np.matmul(U, np.diag(np.sqrt(s))), Vh)\n",
        "    a_tilde = np.matmul(Sigma_sqrt, a)\n",
        "    b_tilde = np.sum(a.flatten() * mu.flatten()) + b\n",
        "    if side == 'positive':\n",
        "        out = norm.cdf(b_tilde / np.linalg.norm(a_tilde))\n",
        "    elif side == 'negative':\n",
        "        out = 1 - norm.cdf(b_tilde / np.linalg.norm(a_tilde))\n",
        "    else:\n",
        "        warnings.warn('Wrong input for argument side.')\n",
        "    \n",
        "    return out"
      ],
      "metadata": {
        "id": "UIJVdLNFrbPM"
      },
      "execution_count": null,
      "outputs": []
    }
  ]
}