"""
Minimal shims for Bio.Application classes removed in Biopython 1.80.

Provides AbstractCommandline, _Option, _Switch, and _Argument with enough
behaviour to support the existing conkit application wrappers.
"""

import subprocess


class _Param:
    def __init__(self, names, description, filename=False, equate=True,
                 is_required=False, checker_function=None):
        self.names = names
        self.description = description
        self.filename = filename
        self.equate = equate
        self.is_required = is_required
        self.checker_function = checker_function
        self.attr = names[1] if len(names) > 1 else names[0]


class _Option(_Param):
    """A -flag value command line option."""
    pass


class _Switch(_Param):
    """A boolean -flag switch (present or absent)."""
    def __init__(self, names, description):
        super().__init__(names, description, equate=False, is_required=False)


class _Argument(_Param):
    """A positional command line argument."""
    def __init__(self, names, description, filename=False, is_required=False,
                 checker_function=None):
        super().__init__(names, description, filename=filename, equate=False,
                         is_required=is_required, checker_function=checker_function)


class AbstractCommandline:
    """Minimal reimplementation of Bio.Application.AbstractCommandline."""

    def __init__(self, cmd, **kwargs):
        self.program_name = cmd
        self._param_map = {p.attr: p for p in self.parameters}
        for attr in self._param_map:
            object.__setattr__(self, attr, None)
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __setattr__(self, name, value):
        if name in ('program_name', 'parameters', '_param_map'):
            object.__setattr__(self, name, value)
            return
        if hasattr(self, '_param_map') and name in self._param_map:
            param = self._param_map[name]
            if param.checker_function is not None and value is not None:
                param.checker_function(value)
        object.__setattr__(self, name, value)

    def _build_args(self):
        args = [self.program_name]
        positional = []
        for param in self.parameters:
            value = getattr(self, param.attr, None)
            if isinstance(param, _Switch):
                if value:
                    args.append(param.names[0])
            elif isinstance(param, _Argument):
                if value is not None:
                    positional.append(str(value))
                elif param.is_required:
                    raise ValueError(f"Required argument {param.attr!r} is not set.")
            else:
                if value is not None:
                    if param.equate:
                        args.append(f"{param.names[0]}={value}")
                    else:
                        args.extend([param.names[0], str(value)])
                elif param.is_required:
                    raise ValueError(f"Required option {param.attr!r} is not set.")
        args.extend(positional)
        return args

    def __str__(self):
        return " ".join(self._build_args())

    def __call__(self):
        result = subprocess.run(
            self._build_args(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"{self.program_name} exited with code {result.returncode}.\n"
                f"stderr: {result.stderr.strip()}"
            )
        return result.stdout, result.stderr
