# Stop a Pipeline

Insure that the you are [prepared](../preparation).

Now all you have to do is run `<your repo root>/tools/pp_stop.sh`.

## pp_stop.sh

Stop the private pipeline started with `pp_start.sh`

It has one optional argument:
1. `<environment profile>` is the set of environment variables that define one pipeline from another. If not provided, it defaults to the username via `${USER}`.

Test for the reader:

1. `pp_stop.sh` uses what `<environment profile>`? [`<repository root>/envs/${USER}`]
1. `pp_stop.sh apple.banana` uses what `<environment profile>`? [`<repository root>/envs/apple.banana`]
1. `pp_stop.sh /cherry/apple/banana` uses what `<environment profile>`? [`/cherry/apple/banana`]


Hence, the most common way to start a private pipeline becomes:

`pp_stop.sh`
