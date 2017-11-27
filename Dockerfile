FROM debian

RUN apt-get update

RUN apt-get install -y \
 samtools=1.3.1-3 \
 python-pip \
 wget

RUN pip install python-Levenshtein numpy scipy

RUN mkdir -p /src; \
 cd /src; \
 wget "https://www.dropbox.com/s/tydhm3ypa7jmbss/scour.py?dl=1" -O scour.py 

RUN mkdir -p /src; \
 cd /src; \
 wget "https://www.dropbox.com/s/o4t5ss6wa37sxrq/split.sh?dl=1" -O split.sh; \
 chmod +x /src/split.sh



