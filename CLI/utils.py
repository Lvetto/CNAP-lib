
def parse_input(user_input):
    t = user_input.split(" ")
    t = [i for i in t if i != " "]
    return t[0], t[1:]

def compile_names(commands):
    return [[i.name] + i.aliases for i in commands]

def find_from_name(command, commands, aliases):
    for n, names in enumerate(aliases):
        if (command in names):
            return commands[n]

    return False

def main_loop(commands, command_names):
    print("Welcome to the interactive analysis tool. Type 'help' for commands.")
    state = {}      # a dict to represent shared memory between commands. Needs to be passed to every command when it's executed
    while True:
        try:
            user_input = input("> ")  # Get user input
            command, args = parse_input(user_input)

            if command == "exit":
                print("Exiting...")
                break
            else:
                command = find_from_name(command, commands, command_names)
                if command:
                    #try:
                        command.execute(state, args)
                    #except:
                    #    print("Something went wrong")
                else:
                    print("Unknown command. Type 'help' for available commands.")

        except (KeyboardInterrupt, EOFError):
            print("\nExiting...")
            break
