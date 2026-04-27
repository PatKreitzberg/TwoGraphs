class MyClass:
    def __init__(self, r1, r2):
        # Use super() to set the initial value to bypass the restriction
        super().__setattr__('r1', r1)
        super().__setattr__('r2', r2)

    def __setattr__(self, name, value):
        if name in ['r1','r2'] and hasattr(self, 'r1'):
            raise AttributeError(f"{name} is immutable")
        super().__setattr__(name, value)


c = MyClass(12, 69)
print(c.r1)
