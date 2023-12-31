FROM ubuntu:22.04
WORKDIR /tmp

# Update and add universe repository
RUN apt-get update && apt-get install -y \
  software-properties-common && \
  add-apt-repository -y universe

# Set timezone to Eastern Timezone
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y tzdata && \
  ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime && \
  dpkg-reconfigure --frontend noninteractive tzdata

COPY system-packages.txt .
COPY requirements.txt .

RUN apt-get update && apt-get install -y git nano $(cat system-packages.txt | perl -pe 's/\n/\ /g')
RUN pip3 install -r requirements.txt

# Clean up extra files
RUN apt-get -y autoremove --purge && apt-get -y clean

WORKDIR /workspace
COPY . .
RUN nodeenv --with-npm /usr/local --force
RUN npm ci

CMD "/bin/bash"