
import itertools

def create_generator():
    return itertools.permutations([1,2,3])

a = create_generator()

for i in a:
    print i
