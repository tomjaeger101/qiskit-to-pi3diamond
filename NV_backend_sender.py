# works at the moment just with python 2.7 
import pika, os
import json
import pickle
import pandas as pa
import sys 
import argparse


"""
    TODO: Maybe use later some argparser functionality, depends on how it will be installed
    def dir_path(string):
        if os.path.isdir(string):
            return string
        else:
            raise NotADirectoryError(string)



    parser_general = argparse.ArgumentParser(description='A parser for the NV backend sender', add_help=False) 

    parser_general.add_argument('--data-dir', type=dir_path, required=True, help='the data directory where data.hdf is located')
    parser_general.add_argument('--user-name', type=str, required=True, help='the user name (to identify to whom the result should be returned)')
    parser_general.add_argument('--timestamp', type=int, required=True, help='The timestamp which is required to uniquely identify a certain measurement')

    args = parser_general.parse_args()
"""

customer_name = "test"
hdf = pa.read_hdf("data.hdf")
timestamp = "123123"

# Access the CLODUAMQP_URL environment variable and parse it (fallback to localhost)
url = os.environ.get('CLOUDAMQP_URL', 'amqps://vsglzbtz:3xiyz9iZcs2S5X-B2hlZ68YlG8n1-5uq@barnacle.rmq.cloudamqp.com/vsglzbtz')
params = pika.URLParameters(url)
#credentials = pika.PlainCredentials('guest', 'guest')


#properties = pika.spec.BasicProperties(user_id='Max') # not in use at the moment
connection = pika.BlockingConnection(params)
channel = connection.channel() # start a channel
queue_result = channel.queue_declare(queue=customer_name) # Declare a queue



channel.basic_publish(exchange='', # exchange if multiple queues are used 
                      routing_key=customer_name,  # queue name 
                      #properties=properties,
                      body=json.dumps({"hdf_file": pickle.dumps(hdf), 'timestamp': timestamp}))

connection.close()