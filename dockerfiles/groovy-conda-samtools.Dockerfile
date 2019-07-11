FROM rsuchecki/samtools:1.9_ea5d3c82fb85c174dce08a4c736ed44c1e6bb7eb

RUN apt-get update \
  && apt-get install -y default-jdk unzip zip \
  && curl -s get.sdkman.io | bash \
  && source "$HOME/.sdkman/bin/sdkman-init.sh" \
  && sdk install groovy