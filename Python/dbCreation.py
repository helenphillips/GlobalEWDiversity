import sqlite3
import sys
import csv

con = sqlite3.connect('C:/Users/hp39wasi/sWorm/EarthwormDatabase/earthwormdb.db')

with con: ## 'with' arguments means all changes are commited    
    cur = con.cursor()    
    cur.execute("CREATE TABLE Authors(AuthorID INTEGER PRIMARY KEY, Title TEXT, Surname TEXT, FirstName TEXT, Email TEXT)");
## If not in this block, need to commit afterwards

## one way of adding rows
cur.execute("INSERT INTO Authors(Title,Surname,FirstName,Email) VALUES('Dr','Doe','Jane','jdoe@ic.ac.uk')");
cur.execute("INSERT INTO Authors(Title,Surname,FirstName,Email) VALUES('Dr','Doe','Jon','jondoe@ic.ac.uk')");
con.commit()

## A safer way of adding rows
title = "Prof."
f_name = "Adam"
s_name = "Smith"
email = "test@nhm.ac.uk"
cur.execute('''INSERT INTO Authors(Title,Surname,FirstName,Email) VALUES(?,?,?,?)''', (title, s_name, f_name, email))
con.commit()

#### Now trying from a CSV file
cur.execute("CREATE TABLE t (Title,Surname,FirstName,Email);") # use your column names here

with open('C:/Users/hp39wasi/sWorm/temp/authortemp.csv','rt') as fin: # `with` statement available in 2.5+
    # csv.DictReader uses first line in file for column headings by default
    dr = csv.DictReader(fin) # comma is default delimiter
    to_db = [(i['Title'], i['Surname'], i['FirstName'], i['Email']) for i in dr]

cur.executemany("INSERT INTO Authors (Title,Surname,FirstName,Email) VALUES (?, ?, ?, ?);", to_db)
con.commit()


### What if CSV contains more than what is needed
with open('C:/Users/hp39wasi/sWorm/temp/authorandmoretemp.csv','rt') as fin: # `with` statement available in 2.5+
    # csv.DictReader uses first line in file for column headings by default
    dr = csv.DictReader(fin) # comma is default delimiter
    test_db = [(i['Title'], i['Surname'], i['FirstName'], i['Email']) for i in dr]

cur.executemany("INSERT INTO Authors (Title,Surname,FirstName,Email) VALUES (?, ?, ?, ?);", test_db)
con.commit()

