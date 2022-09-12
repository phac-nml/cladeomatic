#!/usr/bin/env python3

import sys

tasks = {
    'create': 'Identify population structure and develop typing scheme',
    'benchmark': 'Test developed scheme using labeled samples and scheme',
    'test': 'Test cladeomatic functionality on a small dataset',
    'version': 'Print version and exit',
}

ordered_tasks = [
    'create',
    'benchmark',
]


def print_usage_and_exit():
    print('Usage: cladeomatic <command> [options] <required arguments>', file=sys.stderr)
    print('\nTo get minimal usage for a command use:\ncladeomatic command', file=sys.stderr)
    print('\nTo get full help for a command use one of:\ncladeomatic command -h\ncladeomatic command --help\n', file=sys.stderr)
    print('\nAvailable commands:\n', file=sys.stderr)
    max_task_length = max([len(x) for x in list(tasks.keys())]) + 1
    for task in ordered_tasks:
        print('{{0: <{}}}'.format(max_task_length).format(task), tasks[task], sep=' ', file=sys.stderr)
    sys.exit(0)

def main():

    if len(sys.argv) == 1 or sys.argv[1] in ['-h', '-help', '--help']:
        print_usage_and_exit()

    task = sys.argv.pop(1)

    if task not in tasks:
        print('Task "' + task + '" not recognised. Cannot continue.\n', file=sys.stderr)
        print_usage_and_exit()

    exec('import cladeomatic.' + task)
    exec('cladeomatic.' + task + '.run()')

# call main function
if __name__ == '__main__':
    main()