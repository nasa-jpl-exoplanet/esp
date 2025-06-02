# Development Tools

First, they can be run anywhere. In other words, you can set your PATH to point to this directory and just type `pp_start.sh` to start a private pipeline (pp) or `pp_stop.sh` to stop a private pipeline. You can also not add it to your PATH and just directly access it with a full or relative like `relative/path/to/excalibur/tools/start.sh` or `/absolute/path/to/excalibur/tools/stop.sh`. The scripts will determine where they are and do the right thing.

Second, they can be run on the mentor cluster or remotely on a laptop. It is not recommended to run post2shelve.sh on a laptop as it could take many days to weeks to complete. The major problem with running on a laptop is access to `/proj/sdp/data`. The simplest way to solve the problem is to do these three items:

1. Add a group like excalibur but it can be any name with the GID 1512 to your laptop then add yourself to the group. This step only be done once in the life of your laptop. May also need to log out or reboot your laptop depending on a bunch of factors that are beyond this readme. Just do as much as you need until the name you chose shows up when you execute `groups` at the command line.
1. May need to copy mentor:/etc/ssl/server.* to your laptop /etc/ssl if they do not already exist. Make sure after copy they have the permissions 600 (RW owner only and root should be owner).
1. Use sshfs to mount /proj/sdp/data
   1. Install sshfs if not on your laptop (once in the life of your laptop)
   1. `sudo mkdir -p /proj/sdp/data` (maybe once in the life of your laptop)
   1. `sudo chown ${USER}:${USER}` (maybe once in the life of your laptop)
   1. run sshfs: `sshfs -o allow_other -o default_permissions -o idmap=user ${USER}@excalibur.jpl.nasa.gov:/proj/sdp/data /proj/sdp/data`
   1. `sudo mkdir -p /proj/sdp/${USER}` (maybe once in the life of your laptop)
   1. `sudo chown ${USER}:${USER}` (maybe once in the life of your laptop)
   1. run sshfs: `sshfs -o allow_other -o default_permissions -o idmap=user ${USER}@excalibur.jpl.nasa.gov:/proj/sdp/${USER} /proj/sdp/${USER}`
   
When done, be sure to un-mount `/proj/sdp/data` with `fusermount3 -u /proj/sdp/data`. Same for `/proj/sdp/${USER}`

## Tools

These tools are fairly robust and follow the pythonic guidelines that it should pretty much do out of the box what you want. There may be required values to be set but they will be explicitly noted. Everything else can be left to the default value for the most part. If on your laptop, certainly everything can be left to their defaults.

### autogen.sh

Generate PyXB bindings and update any checksums that go with the new binding. You will normally have to do this when PyTesting fails and it will give you instructions on how to do it.

### inter.sh

At points in time, the operational DB is interred (frozen, copied and moved) in a mausoleum (public facing and accessible data store). The term inter and mausoleum are used partially in humor and partly to describe what is happening to the data. The operational pipeline operates on the concept of current best estimate. To freeze the DB means the data may no longer be the current best estimate; in other words, it has become corpse like in that it resembles the living but is not alive. Once preserved, the data can then be moved to a data store with a more sophisticated UI than the operational system provides; hence, mausoleum.

Interring the data into a mausoleum:

1. list all of the desired values to put into the mausoleum in the file: /proj/sdp/data/nexsci.sv.txt
1. create a directory to place working files and the final result: mkdir -p /proj/sdp/${USER}/R_x.y.z  
    where x, y, and z should be numbers like 1.2.3. Any value of x represents some backward compatibily with all other x of the same value. The value y usually indicates what new features have been added. Lastly, z is patching.
1. run inter.sh to build the mausoleum: tools/inter.sh /proj/sdp/${USER}/R_x.y.z
1. get the file /proj/sdp/${USER}/R_x.y.z.tgz (not the directory) to its caretaker

### make_certs.sh

Generate a pair of TLS certificates to use with EXCALIBUR. The operational pipeline will have to be restarted before your new certs can be used for restricted access functions.

`make_certs.sh /proj/sdp/${USER}/certs`

### post2shelve.sh

Builds a private pipeline database for a private pipeline to use. It requires access to /proj/sdp/data and postgres on the mentor cluster. If using sshfs on a laptop, will also require network access to JPL (VPN).

The personal database will end up in `/proj/sdp/${USER}/db` and named `$USER`.

To run the tool: `post2shelve.sh`

Optionally, `post2shelve.sh` takes one argument to create the database in `/proj/sdp/${USER}/db' and named `$1` (the argument).

### pp_exec.sh

Execute a development task in the private pipeline context started with `pp_start.sh`. Using exec is faster and more resource friendly than using `pp_worker.sh`.

This script has two required and one optional values:
1. `<task name>` is the python package name after excalibur. For instance, if we wanted to run the python packaage `excalibur.cerberus` we would supply `cerberus` as the first argument to the tool.
1. `<target name>` is exactly that. Each star is a target so we would supply the star name as [known by the pipeline](https://excalibur.jpl.nasa.gov:8080/pages/database/targets) like 'GJ 1214'.
1. `<environment profile>` is the set of environment variables that define one pipeline from another. If not provided, it defaults to the username via `${USER}`.

and one optional via an environment variable:
1. `RUNID` is an environment variable that can define the run id for the job you going to execute. If not supplied, the value 17 is used.

Putting it all together:

#### use default run id to run excalibur.eclipse for GJ 1214

`pp_exec.sh eclipse "GJ 1214"`

#### use run id 23 to run excalibur.runtime for all targets

`RUNID=23 pp_exec.sh runtime __all__`

In bash, one can define an environment variable for the life of a command by preceeding the command with the variable definition. In most other shells, `RUNID` would have to be defined separately from executing the command.

### pp_reset.sh

The tool will empty the schedule of jobs, optionally archive the database, and then restart the pipeline. It works with both the operational pipeline and private pipelines..

To reset your private pipeline without archive: `pp_reset.sh`

To reset your private pipeline with archive: `DAWGIE_ARCHIVE=true pp_reset.sh`

To reset operational pipline: `pp_reset.sh ops`

The script triggers off these environment variables:
- `DAWGIE_SFE_PORT` - probably already defined in your bash profile and can be overriden by using the script argument. Default if not defined 8080.
- `EXCALIBUR_HOST` - for private pipelines, use computer name running the pipeline like mentor3. Default if not defined is localhost.
- `DAWGIE_SSL_PEM_MYSELF` - should default to your certificate and will use the one you have defined in your bash profile. Default if not defined /proj/sdp/${EXCALIBUR_USER}/certs/excalibur_identity.pem.

### pp_start.sh

The tool will start a private pipeline. If you are on the mentor machines, you must supply an argument because TCP/IP ports cannot be shared amoung private pipelines. Nothing will force you to do this except an error saying port already in use. On your laptop, can just use the default port.

It has one optional argument:
1. `<environment profile>` is the set of environment variables that define one pipeline from another. If not provided, it defaults to the username via `${USER}`.

Test for the reader:

1. `pp_start.sh` uses what `<environment profile>`? [`<repository root>/envs/${USER}`]
1. `pp_start.sh apple.banana` uses what `<environment profile>`? [`<repository root>/envs/apple.banana`]
1. `pp_start.sh /cherry/apple/banana` uses what `<environment profile>`? [`/cherry/apple/banana`]


Hence, the most common way to start a private pipeline becomes:

`pp_start.sh`

Be aware of these two hidden requirements:
1. DB - the private pipeline database is expected to be in `/proj/sdp/$USER/db` and named `${USER}`. If post2shelve.sh was used to create the DB, then this requirement is already met.
1. runtime - the private pipeline settings and knobs are expected to be the file `/proj/sdp/data/runtime/${USER}.xml`.

### pp_stop.sh

Stop the private pipeline started with `pp_start.sh`

It has one optional argument:
1. `<environment profile>` is the set of environment variables that define one pipeline from another. If not provided, it defaults to the username via `${USER}`.

Test for the reader:

1. `pp_stop.sh` uses what `<environment profile>`? [`<repository root>/envs/${USER}`]
1. `pp_stop.sh apple.banana` uses what `<environment profile>`? [`<repository root>/envs/apple.banana`]
1. `pp_stop.sh /cherry/apple/banana` uses what `<environment profile>`? [`/cherry/apple/banana`]


Hence, the most common way to start a private pipeline becomes:

`pp_stop.sh`

### pp_worker.sh

_Under Construction_

### run.sh

A tool that allows the user to run the github actions on their local computer. The script needs bash and Python 3.12 with pyyaml installed. If on the mentor cluster, use /proj/sdp/.venv/basic.

Two environment variable control the behavior:
- KEEP_CHANGES : when defined, it will keep changes made by black and others
- KEEP_REPORTS : when defined, it will leave the reports in place when successful

You can run the tool anywhere but for example say from the root of the repository like:

```$ tools/run.sh```

To get a list of all of the checks available:

```$ tools/run.sh None```

which will produce output something like this:
```
yaml: dynamic_analysis.yaml
   skipping PyTesting
yaml: static_analysis.yaml
   skipping Style
   skipping Linters
   skipping Security
   skipping Delivery
   skipping Compliance
waiting for jobs to complete
Summary: All test verifications were successful
```

To run just style:

```$ tools/run.sh Style```

To run style and keep the changes that black makes:
```$ KEEP_CHANGES=1 tools/run.sh Style```

Defining the variable as part of the command is a bashism. You can always define them using setenv or export as well but you will lose the easy ability to change your mind later without logging out.

To run more than one check and review the reports:
```$ KEEP_CHANGES=1 KEEP_REPORTS=1 tools/run.sh Style Linters Compliance```


### update_runids.py

Script used by operations to update the information on the UI that indicates what happened in the source code for what run ID. Please just leave this script be.
