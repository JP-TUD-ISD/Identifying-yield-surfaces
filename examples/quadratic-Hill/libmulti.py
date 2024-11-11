import multiprocessing

_Processes = []


def DeleteWorker():
    global _Processes

    # delete old workers
    for worker in _Processes:
        worker.terminate()

    # cleaning
    _Processes = []

    return None


def CreateWorker(NumberOfProcesses):
    '''Create threads and connect queues'''

    # initial work
    global _Processes

    # check for old workers
    if len(_Processes) > 0:
        DeleteWorker()

    # create new workers
    for i in range(NumberOfProcesses):
        _Processes.append(multiprocessing.Process(
                        target=trymp2(i), args=range(10)))

    # start new workers
    for worker in _Processes:
        worker.start()

    return None


def tryingmp():

    print(multiprocessing.current_process())
    # return 1


def trymp2(x):

    print(x)
