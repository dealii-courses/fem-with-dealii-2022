FROM heltai/dealii:vscode

ARG REMOTE_USER=fsi
ARG UID=1000
ARG GID=1000

RUN grep -q $GID /etc/group || groupadd -g $GID $REMOTE_USER
RUN grep -q $UID /etc/passwd || adduser --disabled-password --gecos '' --uid $UID --gid $GID $REMOTE_USER
# make sure we rename the user in case it already existed
RUN [ `id -nu $UID` == "$REMOTE_USER" ] || usermod -l $REMOTE_USER `id -nu $UID`
RUN adduser $REMOTE_USER sudo

#make sure everything is in place
RUN mkdir -p /home/$REMOTE_USER
RUN chown -R $UID:$GID /home/$REMOTE_USER

RUN echo "* hard core 0" >>  /etc/security/limits.conf
RUN echo "* soft core 0" >>  /etc/security/limits.conf
RUN echo "fs.suid_dumpable=0" >> /etc/sysctl.conf
RUN echo "kernel.core_pattern=|/bin/false" >> /etc/sysctl.conf
RUN sysctl -p /etc/sysctl.conf

USER $REMOTE_USER
RUN mkdir -p /home/$REMOTE_USER/.ssh
RUN mkdir -p /home/$REMOTE_USER/.vscode-server-insiders/extensions
RUN mkdir -p /home/$REMOTE_USER/.vscode-server/extensions
ENV HOME /home/$REMOTE_USER
ENV USER $REMOTE_USER
WORKDIR $HOME