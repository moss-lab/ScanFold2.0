FROM ubuntu
#FROM hiroya1024/vienna-rna
#FROM amancevice/pandas
FROM python:3.9.6
COPY . /app
WORKDIR /app
# RUN wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.18.tar.gz
# RUN tar -zxvf ViennaRNA-2.4.18.tar.gz
# WORKDIR ./ViennaRNA-2.4.18
# RUN ./configure --with-python3 --without-kinfold --without-forester --without-kinwalker --without-rnalocmin
# RUN make
# RUN make install
# WORKDIR /app
RUN python3 -m pip install --upgrade pip
RUN pip install -r requirements.txt
#run pip3 install -r requirements.txt
#CMD python ScanFold-Scan2.0.py
ENTRYPOINT ["python", "RunScanFold2.0.py"]
