# the-wizz Python Wrapper

This is a wrapper for the-wizz software (the 2019 version including STOMP maps) by Chris Morrison (see https://github.com/morriscb/the-wizz). 
The code basically wrapps the command line commands into Python classes.


The code is run in the Docker environment provided by the-wizz. 

### Step-by-Step Instructions

#### Download Docker File
First, download the Docker file contining the-wizz. Unless you want to install it yourself. This is not very hard, except the dependencies for STOMP.

```
docker pull morriscb/the-wizz
```

See also instructions here: https://github.com/morriscb/the-wizz


#### Run in Docker

First, create a directory that you share with Docker. This is the directory that can be accessed from the Docker environment. Your door to Docker so to speak.
Let's call this directory the `work` directory:

```
mkdir /path/work
```

where `path` is your favorite path (for example of your project).


Second, enter the Docker environment (not that for this, Docker Desktop has to be running!) by typing:

```
docker run -v /path/work:/home/work/ -it morriscb/the-wizz bash
```

The above line will link the directory `/home/work` **inside** the Docker environment to the above created directory `/path/work`. Whatever you put in that directory can be seen inside the Docker environment and vice versa. **Important is therefore: whatever you want to save to disk from inside the Docker environment, you have to save in _that_ directory**.
The `-it` option lets you run Docker interactively and `bash` starts a bash shell environment.

After executing the above line of code, you are inside the Docker environment. Now, you can run Python as usual. For example, you can run the above code via

```
python WIZZexample.py
```

(obviously after changing the paths and linking files etc.)

To exit the environment, simply type `exit`.

You might need to install new Python packages in your Docker environment. You can do this like you do on your computer, e.g., using `pip`.
For example, the Docker image provided by Chris does not include matplotlib. To install it, simply enter the Docker environment and execute

```
pip install matplotlib
```

Remember that this does **not** change the Docker image. If you exit the environment and restart the it, you will again have to install `matplotlib` (or whatever you installed. **However** there is a trick. You can "re-open" the environment using the following:

First, list all the Docker local "images" on your computer:

```
docker ps -a
```

Then choose the latest one and re-open:

```
docker container start -i <name>
```

where `<name>` is the name of the container.

Basically, this is like resume a `screen` or `byobu` session.


By the way, to check the Docker images, you can type

```
docker images
```


