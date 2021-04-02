array1= [{"test1":21, "test2":22},{"test1":21, "test2":22},{"test1":21, "test2":22},{"test1":21, "test2":22}]
wholelist1=[]
for x in array1:
    templist1=[]
    templist2=[]
    templist3=[]
    templwholelist1=[]
    for values in x.keys():
        templist1+=[values]
    for keys in x.values():
        templist2+=[keys]
    for y in range(len(x)):
        if y <len(x)-1:
            templwholelist1+=[str(templist1[y])+":"+str(templist2[y])+", "]
        if y==len(x)-1:
            templwholelist1+=[str(templist1[y])+":"+str(templist2[y])] 
    wholelist1+=[templwholelist1]

print(wholelist1) 
#print(wholelist1)
