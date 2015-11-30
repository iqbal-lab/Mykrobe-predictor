class Region(object):

    def __init__(self, reference, start, end, forward = True):
        self.reference = reference
        self.start = start
        self.end = end
        self.forward = forward
        if forward:
            self.strand = "forward"
        else:
            self.strand = "backward"

class Gene(Region):

    def __init__(self, name, reference, start, end, forward = True):
        super(self.__class__, self).__init__(reference, start, end, forward)
        self.name = name

