FROM ubuntu:16.04
SHELL ["bash", "-c"]

ENV ILCSOFT /ilcsoft

RUN apt update && apt-get -y install libeigen3-dev build-essential gfortran cmake subversion default-jdk libgsl0-dev libxpm-dev libxft-dev libx11-dev libxext-dev python-dev git wget emacs

RUN mkdir -p $ILCSOFT
RUN git clone -b dev-base https://github.com/eutelescope/ilcinstall $ILCSOFT/ilcinstall
RUN cd $ILCSOFT/ilcinstall && git status

ADD eutelescope.py eutelescope.py

RUN mv eutelescope.py $ILCSOFT/ilcinstall/ilcsoft/eutelescope.py
RUN $ILCSOFT/ilcinstall/ilcsoft-install -i -v $ILCSOFT/ilcinstall/examples/eutelescope/release-standalone.cfg

ENV EUTELESCOPE /ilcsoft/v01-19-02/Eutelescope/master/

WORKDIR ${EUTELESCOPE}

RUN git checkout gbl_work && git pull

RUN chmod +x build_env.sh
#RUN ./build_env.sh
RUN source build_env.sh && cd build && cmake .. && cmake --build . -- -j$(nproc) && make install
#make clean && make install

WORKDIR /ilcsoft/v01-19-02/Eutelescope/master/jobsub/examples/mixed-mode
RUN sed -i 's/\/afs\/desy/\/dummy\/desy/g' duplicated/config.cfg
RUN sed -i 's/\/afs\/desy/\/dummy\/desy/g' standard/config.cfg

RUN mkdir /afs /output /dummy
RUN chmod +rwx /
RUN chmod +rwx /afs
RUN chmod +rwx /output
RUN chmod +rwx /dummy


CMD bash -c " cd /ilcsoft/v01-19-02/Eutelescope/master/ ; source build_env.sh ; cd - ; bash"