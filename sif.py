"""
SIF parser
see http://manual.cytoscape.org/en/stable/Supported_Network_File_Formats.html#sif-format
"""
import re
import io



tab_re = re.compile('.*(\t).*')


class Item:
    def __init__(self, node_a, relation, nodes_b):
        self.node_a = node_a
        self.relation = relation
        self.nodes_b = nodes_b


class SIF:
    def __init__(self, path_or_file):
        # read file and search for tabs
        lines = []

        if not hasattr(path_or_file, 'readlines'):
            fileobj = open(path_or_file)
        else:
            fileobj = path_or_file
        lines, has_tab = self._read_lines(fileobj)
        if has_tab:
            sep = '\t'
        else:
            sep = '\s'
        self.lines = []

        for line in lines:
            items = [x.strip() for x in line.split(sep)]
            node_a = items[0]
            relation = items[1]
            nodes_b = items[2:]
            self.lines.append(Item(node_a, relation, nodes_b))

    def _read_lines(self, fileobj):
        lines = []
        has_tab = False
        for line in fileobj.readlines():
            match = tab_re.match(line)
            if match is not None:
                has_tab = True
            lines.append(line)
        return lines, has_tab

    def __len__(self):
        return len(lines)

