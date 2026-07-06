def opt2cmd(value, cmd_flag):
    """
    Convert a configuration option to a command-line argument string.

    Parameters:
    - value: The value of the configuration option.
    - cmd_flag: The command-line flag corresponding to the option.
    """
    if value is not None:
        return f"{cmd_flag} {value}"
    return ""
