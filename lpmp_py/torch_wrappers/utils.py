import torch
import itertools

def detach(list_of_tensors):
    return [x.detach() for x in list_of_tensors]


def lexico_iter_pairs(lst):
    return itertools.combinations(lst, 2)


def torch_to_numpy_list(list_of_tensors):
    return [x.cpu().detach().numpy() for x in list_of_tensors]


def numpy_to_torch_list(list_of_np_arrays, device, dtype):
    return [torch.from_numpy(x).to(dtype).to(device) for x in list_of_np_arrays]