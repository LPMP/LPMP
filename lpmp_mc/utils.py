import multiprocessing as mp
import sys

def run_func(func, instance_list, parallel = True):
    if parallel:
        try:
            return_dict = run_parallel(func, instance_list)
        except:
            print("Unexpected error: ", sys.exc_info()[0])
            return_dict = run_serial(func, instance_list)
    else:
        return_dict = run_serial(func, instance_list)
    return return_dict

def run_parallel(func, instance_list):
    ctx = mp.get_context('fork')
    manager = ctx.Manager()
    return_dict = manager.dict()
    workers = []
    for b in range(len(instance_list)): 
        worker = ctx.Process(target=func, args=(b, return_dict, instance_list[b]))
        workers.append(worker)
    [w.start() for w in workers]
    for worker in workers:
        worker.join()
        if worker.exitcode != 0:
            raise Exception(f"There was an error during multiprocessing in function {func.__name__} with error code: {worker.exitcode}, in worker: {worker.name}. Possibly too many parallel tasks.")
    return return_dict

def run_serial(func, instance_list):
    return_dict = {}
    for b in range(len(instance_list)): 
        func(b, return_dict, instance_list[b])
    return return_dict