FROM python:2.7
RUN mkdir -p /usr/src/app
WORKDIR /usr/src/app
RUN pip install --upgrade pip
RUN git clone --recursive https://github.com/iqbal-lab/Mykrobe-predictor.git
WORKDIR /usr/src/app/Mykrobe-predictor/mccortex
RUN make    
ENV PATH $PATH:/usr/src/app/Mykrobe-predictor/mccortex/bin
WORKDIR /usr/src/app/Mykrobe-predictor/
RUN pip install mykrobe
CMD mykrobe --help