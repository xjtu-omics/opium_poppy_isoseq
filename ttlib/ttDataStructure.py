class ttUFS(object):
    def __init__(self, setSize):
        self.setSize = setSize
        self.fa = [] 
        self.si = [] 
        for i in range(setSize):
            self.fa.append(i)
            self.si.append(1)

    def find_set(self, x):
        stPoint = x
        while x!=self.fa[x]:
            x = self.fa[x]
        endFa = x
        x = stPoint
        while x!=self.fa[x]:
            tmpX = x
            x =self.fa[x]
            self.fa[tmpX] = endFa
        return endFa
    def union(self, a, b):
        a_father = self.find_set(a) 
        b_father = self.find_set(b)
        if(a_father != b_father):
            a_size = self.si[a_father] 
            b_size = self.si[b_father]
            if(a_size >= b_size): 
                self.fa[b_father] = a_father
                self.si[a_father] = a_size + b_size
                self.si[b_father] = 0
            else:
                self.fa[a_father] = b_father
                self.si[b_father] = a_size + b_size
                self.si[a_father] = 0
