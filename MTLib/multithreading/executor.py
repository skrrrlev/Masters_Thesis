
from threading import Thread, Lock, Event
from time import sleep, time
from multiprocessing import cpu_count

import numpy as np



class Executor:
    '''Executor'''

    def __init__(self, method, *args):

        self.threads: list[Thread] = []
        self.method = method
        self.args = args
        self.janitor_thread: Thread = None
        self.janitor_lock: Lock = None
        self.janitor_event: Event = None
        self.max_threads = max(2,int(cpu_count()*2-3))

    def __enter__(self):
        self.janitor_lock = Lock()
        self.janitor_event = Event()
        self.janitor_thread = Thread(target=self.janitor)
        self.janitor_thread.start()
        
        for current_args in zip(*self.args):
            while True:
                self.janitor_lock.acquire()
                number_of_threads = len(self.threads)
                if number_of_threads<self.max_threads:
                    thread = Thread(target=self.method,args=current_args)
                    thread.start()
                    self.threads.append(thread)
                    self.janitor_lock.release()
                    break
                self.janitor_lock.release()
        
        #print('All threading targets has been started.')

    def janitor(self):
        max_threads = max(2,int(cpu_count()*2-3))
        print(f"Max threads: {max_threads}")

        while True:
            start = time()
            threads: list[Thread] = []
            
            self.janitor_lock.acquire()
            for thread in self.threads:
                if thread.is_alive():
                    threads.append(thread)
            self.threads = threads
            #print(f'Number of active threads: {len(self.threads)}')
            self.janitor_lock.release()
            

            if self.janitor_event.is_set():
                break
            duration = time()-start
            sleep(max(0,0.1-duration))

    def __exit__(self, exc_type, exc_value, traceback):
        
        self.janitor_event.set()
        self.janitor_thread.join()

        for thread in self.threads:
            thread.join()

        
if __name__ == '__main__':
    def method(i,j):
        print(f'Starting thread {i}:{j}.')
        sleep(int(np.random.randint(10)))
        print(f'Ending thread {i}:{j}.')

    with Executor(method,np.arange(1,21),np.arange(1,21)):
        print('Executing')