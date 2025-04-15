# Start a Pipeline

Insure that the you are [prepared](../preparation).

Now all you have to do is run `<your repo root>/tools/pp_start.sh`.

## pp_start.sh

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
