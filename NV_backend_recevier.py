import pika, os
import json 
import sys 
import pickle
import time
import logging

#from Compiler_mcas_qiskit.compiler_qiskit_mcas import construct_full_mcas_file # TODO: change if location is set

logging.basicConfig(filename=os.path.join(os.path.dirname(os.path.realpath(__file__)),'NV_backend_receiver.log'), filemode='w', level=logging.INFO)


output_path = '/Users/maxkeller/Documents/Uni/Hiwi/Wrachtrup_experience/Compiler_and_gates/received_files' # TODO: change if location is set


# Access the CLODUAMQP_URL environment variable and parse it (fallback to localhost)
url = os.environ.get('CLOUDAMQP_URL', 'https://vsglzbtz:3xiyz9iZcs2S5X-B2hlZ68YlG8n1-5uq@barnacle.rmq.cloudamqp.com/vsglzbtz')
params = pika.URLParameters(url)
connection = pika.BlockingConnection(params)
channel = connection.channel() # start a channel 
channel.queue_declare(queue='circuits') # Declare a queue
def callback(ch, method, properties, body):
  json_body = json.loads(body)
  logging.info("qc from {}".format(json_body['username']))
  qc = pickle.loads(bytes.fromhex(json_body['circuit_data'])) # deserialise 
  construct_full_mcas_file(qc, json_body['username'], output_path)
  time.sleep(2) # necessary that no files will be overwritten
  
channel.basic_consume('circuits',
                      callback,
                      auto_ack=True)

logging.info('[*] Waiting for messages:')
channel.start_consuming()
connection.close()