import pubmed_parser
import requests

# requests is a versatile HTTP library in python. Installation: pip install requests

url = "https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/pubmed21n0002.xml.gz"
# url of the file to be downloaded is named as 'url'

response = requests.get(url) # create HTTP response object
# send a HTTP request to the server and save the HTTP response in a response object called 'response'

with open("/Users/sunnie/Desktop/School/Intern/data/pubmed21n0001.xml.gz",'wb') as f:
# Saving received content as a xml.gz file in binary format
#   
    f.write(response.content)
# write the contents of the response(response.content) to a new file in binary mode.

artisle_list = pubmed_parser.parse_medline_xml("/Users/sunnie/Desktop/School/Intern/data/pubmed21n0001.xml");
# https://titipata.github.io/pubmed_parser/api.html parse MEDLINE XML file
# .xml or .xml.gz doesn't matter

print(artisle_list[0].get('title'));
