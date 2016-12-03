FROM python:2.7
RUN mkdir -p /usr/src/app
WORKDIR /usr/src/app
RUN pip install --upgrade pip
COPY . /usr/src/app
RUN pip install -r requirements.txt 
RUN python setup.py install
CMD mykrobe --help