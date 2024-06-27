import pymysql
from pymysql.constants import CLIENT

conn = pymysql.connect(
    host="localhost",
    user="root", password="",
    database="crisper_guide",
    charset="utf8",
    client_flag=CLIENT.MULTI_STATEMENTS)