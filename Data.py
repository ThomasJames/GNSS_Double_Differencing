class Node(object):
    def __init__(self, d, n=None):
        self.data = d
        self.next_node = n

    def get_next(self):
        return self.next_node

    def set_next(self, n):
        self.next_node = n

    def get_data(self):
        return self.data

    def set_data(self, d):
        self.data = d

class LinkedList(object):
    def __init__(self, r=None):
        self.root = r
        self.size = 0

    def get_size(self):
        return self.size

    def add(self, d):
        new_node = Node(d, self.root)
        self.root = new_node
        self.size += 1


# HOW TO REMOVE A NODE
    def remove(self, d):
        this_node = self.root
        prev_node = None

        while this_node:
            if this_node.get_data() == d:
                if prev_node:
                    prev_node.set_next(this_node.get_next())
                else:
                    self.root = this_node
                self.size -= 1
                return True # REMOVE THE DATA
            else:
                prev_node = this_node
                this_node = this_node.get_next()
        return False # DATA NOT FOUND


    def find(self, d):
        this_node = self.root
        while this_node:
            if this_node.get_data() == d:
                return d
            else:
                this_node = this_node.get_next()
        return None

    
my_list = LinkedList()
my_list.add(6)
my_list.add(7)
my_list.add(10)

print(my_list.find(6))

                
        

      

      


    

    








# def romanToInt(s):
#     l = []
#     for i in s:
#         l.append(i)
#     print(l)
#
#     real = []
#     for i in range(len(l)):
#         if l[i] == "M":
#             real.append(1000)
#         elif l[i] == "D":
#             real.append(500)
#         elif l[i] == "C":
#             real.append(100)
#         elif l[i] == "L":
#             real.append(50)
#         elif l[i] == "X":
#             real.append(10)
#         elif l[i] == "V":
#             real.append(5)
#         elif l[i] == "I":
#             real.append(1)
#
#
#     real = sum(real)
#     counter = 0
#
#     for i in range(len(l)):
#         try:
#             if l[i] == "I" and l[i+1] == "V":
#
#
#         except:
#             continue
#         try:
#             if l[i] == "I" and l[i + 1] == "X":
#                 counter += 1
#         except:
#             continue
#         try:
#             if l[i] == "X" and l[i + 1] == "L":
#                 counter += 40
#         except:
#             continue
#         try:
#             if l[i] == "X" and l[i + 1] == "C":
#                 counter += 40
#         except:
#             continue
#         try:
#             if l[i] == "C" and l[i + 1] == "D":
#                 counter += 100
#         except:
#             continue
#         try:
#             if l[i] == "C" and l[i + 1] == "M":
#                 counter += 100
#         except:
#             continue
#
#     print(counter)
#
#     return real - counter
#
# a = "IX"
# print(romanToInt(a))
#