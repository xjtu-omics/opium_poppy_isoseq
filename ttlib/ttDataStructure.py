class ttUFS(object):
    # 初始化
    def __init__(self, setSize):
        # 初始化两个列表
        self.setSize = setSize
        self.fa = [] # 保存元素所属集合的代表元素
        self.si = [] # 保存父节点包含的元素个数
        for i in range(setSize):
            self.fa.append(i)
            self.si.append(1)

    # FIND-SET操作
    # 采用递归的策略定位父节点
    # 在父节点查找过程中，将当前节点连接到父节点上，进行路径压缩
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

    # UNION操作
    # 将a和b两个集合合并在一起
    def union(self, a, b):
        a_father = self.find_set(a) # 获取两元素所在集合的代表元素
        b_father = self.find_set(b)
        if(a_father != b_father):
            a_size = self.si[a_father] # 获取两元素所在集合的大小
            b_size = self.si[b_father]
            if(a_size >= b_size): # 将规模较小的集合合并到规模较大的集合下面
                self.fa[b_father] = a_father
                self.si[a_father] = a_size + b_size
                self.si[b_father] = 0
            else:
                self.fa[a_father] = b_father
                self.si[b_father] = a_size + b_size
                self.si[a_father] = 0
