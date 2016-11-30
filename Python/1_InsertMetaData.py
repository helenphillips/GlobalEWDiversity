import sqlite3
import sys
import csv

with open('C:/Users/hp39wasi/sWorm/EarthwormDatabase/Singh_2016/Singh_2016 - MetaData.csv','rt') as fin: # `with` statement available in 2.5+
    # csv.DictReader uses first line in file for column headings by default
    dr = csv.DictReader(fin) # comma is default delimiter
    title = [(i['Article_Title']) for i in dr] ## That works...but is it how I want to do it
    

    test_db = [(i['Article_Title'], i['Article_Year'], i['Article_FirstAuthorSurname'], i['Article_Journal'],i['Article_DOI'], i['DataProvider_Title'], i['DataProvider_Surname'], i['DataProvider_FirstName'], i['DataProvider_MiddleInitials'], i['DataProvider_Email'], i['DataProvider_Institute'], i['DataProvider_Department'], i['Number_of_Studies'], i['Total_Number_ofSites'], i['Total_Number_ofSpecies'], i['Notes'], i['BibKey']) for i in dr]





cur.execute("SELECT AuthorID FROM Authors WHERE Surname = ? AND FirstName = ?", (s_name2, f_name))
data=cur.fetchone()
   
if data is None:
    print('There is no person named %s %s, inserting'%(f_name,s_name2))
    cur.execute('''INSERT INTO Authors(Title,Surname,FirstName,Email) VALUES(?,?,?,?)''', (title, s_name2, f_name, email))
    con.commit()
else:
    print('%s %s found with rowid %s'%(f_name,s_name2,data[0]))

																
