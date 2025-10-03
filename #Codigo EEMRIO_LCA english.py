#Codigo EEMRIO by pilinya

#This information is necessary in Python, for activate the funtions 
#If they are not installed is necessary to do it, write pip install XXXXX, being in the floor where are pip.py (in my case C:\Users\LIP\AppData\Local\Programs\Python\Python311\Scripts)

#Activate matrix
import pandas as pd #To real CSV and dataframes
import numpy as np #To work with matrix
from scipy import io, integrate, linalg, signal 
from scipy.sparse.linalg import cg, eigs 


#To create list
from typing import List, Tuple, Dict, Callable, Iterable, Union

#To paint maps
import plotly.express as px

#This is to read input-output matrix, you have to be in the folder where CSVs are

df_InOut = pd.read_csv('io.csv', sep=';') #input output matriz
df_exten= pd.read_csv('extension.csv', sep=';') #Total emissions for each activity
df_exten_social=pd.read_csv('extension_social.csv', sep=';')#Social information for each activity
df_X = pd.read_csv('X.csv', sep=';') #Vector X, total units of each activity

InOut=df_InOut.to_numpy() #To transform the dataframes (read_csv stores as dataframe) to matrix
exten=df_exten.to_numpy()
exten_social=df_exten_social.to_numpy()
X=df_X.to_numpy()

#To calculate the unitary impacts of each activity (matrix R)

(a,b)=np.shape(exten)
R=np.zeros((a-2,b)) #create a matrix of zeros to store the new values
R=np.vstack([exten[:2,:],R]) #To put the headings (the rows)

val_exten=np.delete(exten,(0,1),axis=0) #To remove the headings and be able to divide matrix
val_exten=np.delete(val_exten,(0,1,2),axis=1)
val_exten=val_exten.astype(float) #Numbers are stored as string, to change it to numbers
val_X=np.delete(X,(0,1),axis=0) #Same as before
val_X=np.delete(val_X,(0,1,2),axis=1)
val_X=val_X.astype(float)

R_vals=val_exten/val_X #Calculate the value of R (total emissions/total units)
R_vals=np.nan_to_num(R_vals) #Remove the NAN values (some values are 0/0)
R[2:,3:]=R_vals

R[:,:3]=exten[:,:3] #To put the headings (the columns)

#To do R, but with social information
(a1,b1)=np.shape(exten_social)
R_soc=np.zeros((a1-2,b1)) #create a matrix of zeros to store the new values
R_soc=np.vstack([exten_social[:2,:],R_soc]) #To put the headings (the rows)

val_exten_s=np.delete(exten_social,(0,1),axis=0) #To remove the headings and be able to divide matrix
val_exten_s=np.delete(val_exten_s,(0,1,2),axis=1)
val_exten_s=val_exten_s.astype(float) #Numbers are stored as string, to change it to numbers


R_vals_soc=val_exten_s/val_X #Calculate the value of R (total social impact/total units)
R_vals_soc=np.nan_to_num(R_vals_soc) #Remove the NAN values (some values are 0/0)
R_soc[2:,3:]=R_vals_soc

R_soc[:,:3]=exten_social[:,:3] #To put the headings (the columns)

#To calculate the M matrix (impact value)
country= input ('Which country or region is?(In capital letters based on the code ISO 3166-1 alpha2): ') #To select data
activity=input ('Which activity is? (Exactly the same name of EXIOBASE, with same capital and small letters): ')
df_countries = pd.read_csv('countries.csv', sep=';') #To read the countries
df_activities=pd.read_csv('activities.csv', sep=';') #To read the activities

countries_list=list(df_countries) #To transform the dataframes (read_csv stores as dataframe) to a list
activities_list=list(df_activities)

coun=countries_list.index(country) #To look where is the country you look for
act=activities_list.index(activity) #To look where is the activity you look for
index=3+coun*164+act#The index of the column of Input-output matrix

equis=InOut[2:,index] #To get the numbers of IO with outheadings
equis=equis.astype(float) #Numbers are stored as string, to change it to numbers

R_val=R[2:,3:] #To get the numbers of R without headings
R_val_soc=R_soc[2:,3:] #To get the numbers of social R without headings

M=R_val*equis.T #Calculate the vector M (impact vector), element-by-element multiplication to get impacts by country

#To calculate equivalent CO2

CF=[1, 298, 36.8, 1507, 1281.17, 26100] #This factors are based on EF 3.0, change if you want another, change it (CO2, N2O, CH4, HFC, PFC, SF6)
M_eq_val=CF@M

M_eq=InOut[2:,:3].T #I do the transposed matrix, because I don't know how to add columns, I only know how to add rows 
M_eq=np.vstack([M_eq,M_eq_val]) #To add headings in M_eq

np.savetxt('M_eq1ZAm1.txt',M_eq, fmt="%s") #To save the M_eq as txt

#Use of ecoinvent to update and calculate and specific impact

P_eq=list(pd.read_csv('countries.csv', sep=';')) #Headings of P_eq
activities=list(pd.read_csv('activities.csv', sep=';')) #To obtain the activities

M_eq_val_sum=np.zeros(len(P_eq)) #Create a matrix of zeros to store values

for i in range (len(P_eq)):
    M_eq_val_sum[i]=sum(M_eq_val[i*164:i*164+len(activities)]) #Total impacts per country


imp_total=sum(M_eq_val_sum)#To calculate the percentage of impact per country
M_eq_pais_porc=M_eq_val_sum/abs(imp_total)

#This three lines only if you have information, unless remove
#imp_ese_pais=input('Which is the unitary impact for this contry (decimals with point): ') 
#imp_global=input('Which is the unitary impact globally (decimals with point): ')
#param_k=imp_ese_pais/abs(imp_global) #k value

valor_LCA_conv= input ('Which is the update and specific impact of this activity? Global or regional obtained from ecoinvent or other database (decimals with point): ') #To obtain an update and specific impact
valor_LCA_conv=float(valor_LCA_conv) #Numbers are stored as string, to change it to numbers

P_eq_val=M_eq_pais_porc*valor_LCA_conv  #To calculate updated and specific impact per country
P_eq=np.vstack([P_eq,P_eq_val]) #To add headings in P_eq (The rows)

#Social impact
M_social_val=R_val_soc*equis.T #Calculate the vector M (social impact vector), element-by-element multiplication to get impacts by country
M_social=list(pd.read_csv('countries.csv', sep=';')) #Headings of M_social
activities=list(pd.read_csv('activities.csv', sep=';')) #To obtain the activities

M_social_2=np.zeros((np.size(M_social_val,0),len(M_social))) #Create a matrix of zeros to store values

for i in range (len(M_social)):
    M_social_2[0][i]=sum(M_social_val[0][i*164:i*164+len(activities)]) #Sum per country social impact 1
    M_social_2[1][i]=sum(M_social_val[1][i*164:i*164+len(activities)]) #Sum per country social impact 2
    M_social_2[2][i]=sum(M_social_val[2][i*164:i*164+len(activities)]) #Sum per country social impact 3
    M_social_2[3][i]=sum(M_social_val[3][i*164:i*164+len(activities)]) #Sum per country social impact 4
    M_social_2[4][i]=sum(M_social_val[4][i*164:i*164+len(activities)]) #Sum per country social impact 5
    M_social_2[5][i]=sum(M_social_val[5][i*164:i*164+len(activities)]) #Sum per country social impact 6

value_unit=val_X[0][index-3]
M_social_2=np.divide(M_social_2,value_unit)
M_social=np.vstack([M_social,M_social_2]) 

#To save the values as TXT (change the name of file if you want)
np.savetxt('P_eq_AU_COBRE.txt',P_eq, fmt="%s") #To save P_eq
np.savetxt('M_social_pol_cobre.txt',M_social, fmt="%s") #To save the M_social

#To do the maps of impacts (I don't why these lines sometimes do not work, re-run the code if maps are not painted). Change P_eq per M_social to paint social values

countries_alpha3=pd.read_csv('countries_alpha3.csv', sep=';') #To create the dataframe
P_eq=pd.DataFrame(P_eq.T)
P_eq=P_eq.set_axis(['EXIOBASE','Impact'],axis=1)

P_eq_todo=pd.merge(countries_alpha3,P_eq,how="outer") #Join two dataframe in one
P_eq_todo=P_eq_todo.dropna(axis=1) #To eliminate empty columns (something should not work correctly)
P_eq_todo['Impact']=P_eq_todo['Impact'].astype(float)
#En caso de querer guardarlo: np.savetxt('P_eq_todo_ZAm1.txt',P_eq_todo, fmt="%s") 

fig= px.choropleth(P_eq_todo,locations='ISO3',color='Impact',color_continuous_scale=px.colors.sequential.Plasma) #To paint the map, copy from internet


fig.show() #To show the map
