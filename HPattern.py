#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import imageio


# In[2]:


def MakeCircle(VH,vh,rad):
    #VH are the size of the image.
    #vh is the location of the center
    #rad is the radius of the circle
    a,b = np.indices(VH)
    v,h = vh # location of center
    a -= v
    b -= h
    dist = np.sqrt(a*a + b*b)
    answ = (dist<rad).astype(int)
    return answ


# In[3]:


def HbondPaternMG(ltrs): #ltrs is the dna string
    l= len(ltrs)
    mg = np.zeros((125,1,3),int)
    for i in range(l):
        ans = np.zeros((125,25,3), int)
        if ltrs[i] == 'c':
            # the pattern of 'c'
            mask1 = MakeCircle((125,25),(50,12),5 )
            mask2 = MakeCircle((125,25),(75,12),5 )
            mask3 = MakeCircle((125,25),(100,12),5 )
            ans[:,:,2] = mask1 * 255
            ans[:,:,0] += mask2 * 255
            ans[:,:,0] += mask3 *255
            ans[:,24,:]= 255
        if ltrs[i]== 'g':
            # the pattern of 'g'
            mask1 = MakeCircle((125,25),(25,12),5 )
            mask2 = MakeCircle((125,25),(50,12),5 )
            mask3 = MakeCircle((125,25),(75,12),5 )
            ans[:,:,0] = mask1 * 255
            ans[:,:,0] += mask2 * 255
            ans[:,:,2] = mask3 *255
            ans[:,24,:]= 255
        if ltrs[i]== 't':
            # the pattern of 't'
            mask1 = MakeCircle((125,25),(25,12),5 )
            mask2 = MakeCircle((125,25),(50,12),5 )
            mask3 = MakeCircle((125,25),(75,12),5 )
            mask4 = MakeCircle((125,25),(100,12),5 )
            ans[:,:,0] = mask1 * 255
            ans[:,:,1] = mask1 * 255
            ans[:,:,2] = mask1 * 255
            ans[:,:,0] += mask2 * 255
            ans[:,:,2] += mask3 *255
            ans[:,:,0] += mask4 *255
            ans[:,24,:]= 255
        if ltrs[i]== 'a':
            # the pattern of 'a'
            mask1 = MakeCircle((125,25),(25,12),5 )
            mask2 = MakeCircle((125,25),(50,12),5 )
            mask3 = MakeCircle((125,25),(75,12),5 )
            mask4 = MakeCircle((125,25),(100,12),5 )
            ans[:,:,0] = mask1 * 255
            ans[:,:,2] = mask2 * 255
            ans[:,:,0] += mask3 *255
            ans[:,:,0] += mask4 *255
            ans[:,:,1] = mask4 *255
            ans[:,:,2] += mask4 *255
            ans[:,24,:]= 255
        a1 = ans[:,:,0] ==0
        a2 = ans[:,:,1] ==0
        a3 = ans[:,:,2] == 0
        mask = a1*a1*a3
        # the back ground is grey
        for p in range(3):
            ans[:,:,p] = mask*170 + (1-mask)*ans[:,:,p]
        # collect the patterns into one image
        mg= np.concatenate(( mg,ans),1) 
    return mg


# In[4]:


# to create a matrix for the pattern of the the H-bonds
def HbondPatern(ltrs): #ltrs is the dna string
    l= len(ltrs)
    hmat= np.zeros((4,l),'<U1')
    for i in range(l):
        if ltrs[i] == 'c':
            hmat[:,i]= '-','a','d', 'd' 
        if ltrs[i]== 'g':
            hmat[:,i]= 'd','d','a','-'
        if ltrs[i]== 't':
            hmat[:,i]= 'n','d','a','d'
        if ltrs[i]== 'a':
            hmat[:,i]= 'd','a','d','n'
        if ltrs[i]== '.':
            hmat[:,i]= '.','.','.','.'
    return hmat


# In[5]:


# get the matrix of the match
# to compare 2 sequences together and get the match of the H-bonds
def Compare(ltrs, conseq):
    p1= HbondPatern(ltrs)
    p2= HbondPatern(conseq)
    l= len(ltrs)
    comp= np.zeros((4,l),'<U1')
    for i in range(4):
        for j in range(l):
            if p1[i,j]== p2[i,j]:
                comp[i,j]= p1[i,j]
            else:
                comp[i,j]= ' '
    return comp


# In[6]:


# to get the percetage of the H-bond pattern match between 2 sequences
def PercentMatch(ltrs, conseq):
    p1= HbondPatern(ltrs)
    p2= HbondPatern(conseq)
    l= len(ltrs)
    comp= np.zeros((4,l)) #the comparison matrix
    for i in range(4):
        for j in range(l):
            if p1[i,j]== p2[i,j]:
                comp[i,j]= 1
            else:
                comp[i,j]= 0
    mmatch= comp.sum() #sum of the comp matrix
    compsize= comp.size # number of the elements in the comp matrix
    percent= (mmatch/compsize)*100 # match %
    return percent


# In[7]:


# to get the image of the compared pattern
def HbondMG(mat): # mat is the comparison matrix
    rows,clmns= mat.shape
    mg = np.zeros((125,1,3),int)
    for i in range(clmns):
        for j in range (rows):
            ans = np.zeros((125,25,3), int)
            if mat[j,i]== ' ':
                ans[:,24,:]= 255
            if mat[0,i] == 'd':
                #print('a')
                mask = MakeCircle((125,25),(25,12),5 )
                ans[:,:,0] += mask * 255 
                ans[:,24,:]= 255
            if mat[0,i] == 'a':
                #print('b')
                mask = MakeCircle((125,25),(25,12),5 )
                ans[:,:,2] += mask * 255 
                ans[:,24,:]= 255
            if mat[0,i] == 'n':
                #print('b')
                mask = MakeCircle((125,25),(25,12),5 )
                ans[:,:,0] += mask * 255 
                ans[:,:,1] += mask * 255 
                ans[:,:,2] += mask * 255 
                ans[:,24,:]= 255
            if mat[1,i] == 'd':
                #print('a')
                mask = MakeCircle((125,25),(50,12),5 )
                ans[:,:,0] += mask * 255 
                ans[:,24,:]= 255
            if mat[1,i] == 'a':
                #print('b')
                mask = MakeCircle((125,25),(50,12),5 )
                ans[:,:,2] += mask * 255 
                ans[:,24,:]= 255
            if mat[1,i] == 'n':
                #print('b')
                mask = MakeCircle((125,25),(50,12),5 )
                ans[:,:,0] += mask * 255 
                ans[:,:,1] += mask * 255 
                ans[:,:,2] += mask * 255 
                ans[:,24,:]= 255
            if mat[2,i] == 'd':
                #print('a')
                mask = MakeCircle((125,25),(75,12),5 )
                ans[:,:,0] += mask * 255 
                ans[:,24,:]= 255
            if mat[2,i] == 'a':
                #print('b')
                mask = MakeCircle((125,25),(75,12),5 )
                ans[:,:,2] += mask * 255 
                ans[:,24,:]= 255
            if mat[2,i] == 'n':
                #print('b')
                mask = MakeCircle((125,25),(75,12),5 )
                ans[:,:,0] += mask * 255 
                ans[:,:,1] += mask * 255 
                ans[:,:,2] += mask * 255 
                ans[:,24,:]= 255
            if mat[3,i] == 'd':
                #print('a')
                mask = MakeCircle((125,25),(100,12),5 )
                ans[:,:,0] += mask * 255 
                ans[:,24,:]= 255
            if mat[3,i] == 'a':
                #print('b')
                mask = MakeCircle((125,25),(100,12),5 )
                ans[:,:,2] += mask * 255 
                ans[:,24,:]= 255
            if mat[3,i] == 'n':
                #print('b')
                mask = MakeCircle((125,25),(100,12),5 )
                ans[:,:,0] += mask * 255 
                ans[:,:,1] += mask * 255 
                ans[:,:,2] += mask * 255 
                ans[:,24,:]= 255
            a1 = ans[:,:,0] ==0
            a2 = ans[:,:,1] ==0
            a3 = ans[:,:,2] == 0
            mask = a1*a1*a3
            for p in range(3):
                ans[:,:,p] = mask*170 + (1-mask)*ans[:,:,p]
            #print('e')
        mg= np.concatenate(( mg,ans),1) 
    return mg


# In[8]:


# to compare 3 sequences together and to get the comparison matrix as a result
def Compare2(ltrs, conseq1, conseq2):
    p1= HbondPatern(ltrs)
    p2= HbondPatern(conseq1)
    p3= HbondPatern(conseq2)
    l= len(ltrs)
    comp= np.zeros((4,l),'<U1')
    for i in range(4):
        for j in range(l):
            if p1[i,j]== p2[i,j] and p1[i,j]==p3[i,j]:
                comp[i,j]= p1[i,j]
            else:
                comp[i,j]= ' '
    return comp


# In[9]:


#comparison of 2 matrecies together
def CompareMats(mat1, mat2):
    rows,clmns = mat1.shape
    comp= np.zeros((rows,clmns),'<U1')
    for i in range(rows):
        for j in range(clmns):
            if mat1[i,j]== mat2[i,j]:
                comp[i,j]= mat1[i,j]
            else:
                comp[i,j]= ' '
    return comp


# In[ ]:




