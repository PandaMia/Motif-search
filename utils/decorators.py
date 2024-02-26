from functools import wraps


def error_handler(func):
    @wraps(func)
    def handled_function(*args, **kwargs) -> dict:
        try:
            return func(*args, **kwargs)
        except Exception as e:
            error = f"""LOCATION: Module - {func.__module__}.py, Function - {func.__name__} 
                        ERROR: {e}""".replace("  ", "")
            print(error)
            response = {"error": error}
            return response
    return handled_function