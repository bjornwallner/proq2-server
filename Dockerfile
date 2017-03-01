FROM ubuntu:16.04
RUN apt update && apt -y install git wget gnuplot emboss
RUN git clone http://github.com/bjornwallner/proq2-server/
RUN cd proq2-server/ && git pull
RUN cd proq2-server/ && bin/download_dbs.sh
RUN mkdir /proq2-server/server/output/
VOLUME /proq2-server/server/output/
ENTRYPOINT ["/proq2-server/bin/run_proq2"]
