import sqlite3
import sys
import csv

con = sqlite3.connect('C:/Users/hp39wasi/sWorm/EarthwormDatabase/earthwormdb.db')

## Need an if statement, to check whether table has already been created
with con: ## 'with' arguments means all changes are commited    
    cur = con.cursor()    
    cur.execute("CREATE TABLE IF NOT EXISTS Institutes(instituteID INTEGER PRIMARY KEY, name TEXT, department TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS Dataproviders(dataprovidersID INTEGER PRIMARY KEY, title TEXT, surname TEXT, firstname TEXT, middleinitials TEXT, email TEXT, instituteID INTEGER, FOREIGN KEY(instituteID) REFERENCES Institutes(instituteID))")
    cur.execute("CREATE TABLE IF NOT EXISTS Journals(journalID INTEGER PRIMARY KEY, name TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS Articles(articleID INTEGER PRIMARY KEY, title TEXT, doi TEXT, year DATE, firstauthor_surname TEXT, dataproviderID INTEGER, numberofstudies INTEGER, numberofsites INTEGER, numberofspecies INTEGER, Bibkey TEXT, journalID INTEGER, FOREIGN KEY(dataproviderID) REFERENCES Dataproviders(dataprovidersID), FOREIGN KEY(journalID) REFERENCES Journals(journalID))");

