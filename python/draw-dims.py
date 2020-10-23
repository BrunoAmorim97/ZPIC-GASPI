dims = [4,4]

proc = 0

for row in range(dims[0]):
    for column in range(dims[1]):

        print(" {:>2} ".format(proc), end="")
        proc += 1

        if(column != dims[0] - 1):
            print("|", end="")
        
    print()
    
    if(row != dims[1] - 1):
        print("-----" * dims[0])
