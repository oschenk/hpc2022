# Introduction to the ICS cluster

Overview of the cluster https://intranet.ics.usi.ch/HPC
* hardware and software configuration
* scientific software
* access the infrastructure to run jobs

**Additional reading:**
* Remote machines section of the MIT lecture https://missing.csail.mit.edu/2020/command-line/


## Connect to the cluster

You have received the username and password to your USI email. With these
you can connect to the cluster using the command below.
If you have not received the credentials, please contact the TA.

```
$ ssh studXX@hpc.ics.usi.ch
```

Use command `exit` or `ctrl+D` to disconnect from the session.

## Set-up Public Key Authentication for SSH

### Generate ssh-key

To avoid typing the password every time you can generate a ssh-key on your laptop and copy it to the cluster.

Run on your laptop:

```
$ ssh-keygen
$ ssh-copy-id -i ~/.ssh/id_rsa.pub studXX@hpc.ics.usi.ch
```

### Add host configuration to `~/.ssh/config`

Add the following configuration to the file `~/.ssh/config` on your laptop. If the file does not exist, create it.

```
Host icsmaster
   Hostname hpc.ics.usi.ch
   Port 22
   User studXX
   IdentityFile ~/.ssh/id_rsa
```

### Done!
Now you can log-in with the ssh-key instead of password.

```
$ ssh icsmaster
```
