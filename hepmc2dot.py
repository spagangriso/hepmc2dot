#!/usr/bin/env python

import math
import re # regex
import numpy as np


def _get_dot_particle(prod_vtx_barcode, end_vtx_barcode,
                      particle_barcode, particle_id, particle_energy):
    """
    Returns a string containing a DOT formatted edge which represents a particle travelling from
    the given production to the given end vertex. If end_vtx_barcode is None, the edge will connect
    to a dummy end vertex.
    """
    prod_vtx = _get_node_name(prod_vtx_barcode)
    if not end_vtx_barcode:
        end_vtx = _get_node_name(particle_barcode, is_dummy=True)
    else:
        end_vtx = _get_node_name(end_vtx_barcode)

    particle_dot = '    {prod_vtx} -> {end_vtx} [label="p #{bc}\\n' \
                   'id={part_id}\\n' \
                   'E={energy:.0f}"];\n'.format(prod_vtx=prod_vtx,
                                                end_vtx=end_vtx,
                                                bc=particle_barcode,
                                                part_id=particle_id,
                                                energy=float(particle_energy))
    return particle_dot


def _get_node_name(barcode, is_dummy=False):
    """
    Returns the DOT node name for the given barcode

    DOT nodes represent either an interaction vertex or a dummy end vertex for a final state
    particle
    """
    if is_dummy:
        dummy = 'dummy_'
    else:
        dummy = ''

    abs_bc = abs(int(barcode))
    node_name = 'V_{dummy}{abs_bc}'.format(dummy=dummy, abs_bc=abs_bc)
    return node_name


def _get_dot_vertex(barcode, r, z, is_dummy=False, scale=1.):
    """
    Generates a DOT formatted string representing an interaction (or dummy) vertex for the given
    barcode and coordinates
    """
    barcode = int(barcode)
    r = float(r)
    z = float(z)
    if is_dummy:
        attrib = 'shape=none,label=""'
    else:
        attrib = r'label="vtx #{bc}\nr={r:.2f},z={z:.2f}"'.format(bc=barcode, r=r, z=z)

    vtx_name = _get_node_name(barcode, is_dummy)
    dot = '    {node} [{attrib},pos="{zpos:.3f},{rpos:.3f}!"];\n'.format(node=vtx_name,
                                                                         attrib=attrib,
                                                                         zpos=z*scale,
                                                                         rpos=r*scale)
    return dot


class HepDotWriter(object):
    """
    Generates a dot file representing the given particles, interaction vertices and events
    """

    def __init__(self, dotfile): #TODO: implement next:, vtx_threshold=np.nan, scale=1.):
        self.dotfile = open(dotfile, 'w')

        self.event_open = False
        self.cur_vtx_barcode = None
        self.cur_vtx_r = None
        self.cur_vtx_z = None

        # TODO: implement next:
        #self.vtx_threshold = np.nan # vtx_threshold
        self.scale = 1. #scale
        # primary only:
        #self.vtx_threshold = 200000
        #self.scale = 50.
        # all particles:
        #self.vtx_threshold = np.nan
        #self.scale = 2.

    def newEvent(self, raw_hepmc_line):
        self._end_opened_event()
        self._begin_event(raw_hepmc_line)

    def newVertex(self, raw_hepmc_line):
        hepmc = raw_hepmc_line.split()
        vtx_barcode_column = 1
        vtx_barcode = int(hepmc[vtx_barcode_column])
        vtx_abs_barcode = abs(vtx_barcode)

        # TODO: implement next
        #if vtx_abs_barcode > self.vtx_threshold:
        #    return

        x_column = 3
        x = float(hepmc[x_column])
        y_column = 4
        y = float(hepmc[y_column])
        z_column = 5

        self.cur_vtx_z = float(hepmc[z_column])
        self.cur_vtx_r = math.sqrt(x**2 + y**2)

        self.cur_vtx_barcode = vtx_barcode
        dot_vtx = _get_dot_vertex(vtx_barcode,
                                  self.cur_vtx_r,
                                  self.cur_vtx_z,
                                  scale=self.scale)
        self.dotfile.write(dot_vtx)

    def addOutgoingParticle(self, raw_hepmc_line):
        line = raw_hepmc_line.split()

        particle_barcode_column = 1
        particle_barcode = line[particle_barcode_column]

        particle_id_column = 2
        particle_id = line[particle_id_column]

        particle_energy_column = 6
        particle_energy = float(line[particle_energy_column])

        mom_x_column = 3
        mom_x = float(line[mom_x_column])
        mom_y_column = 4
        mom_y = float(line[mom_y_column])
        mom_z_column = 5
        mom_z = float(line[mom_z_column])

        end_vtx_barcode_column = 11
        end_vtx_barcode = abs(int(line[end_vtx_barcode_column]))

        # TODO: implement next
        #if end_vtx_barcode > self.vtx_threshold:
        #    return

        if not end_vtx_barcode:
            # create dummy end node for partiles that don't have end vertices

            mom_r = math.sqrt(mom_x**2 + mom_y**2)
            mom_abs = math.sqrt(mom_r**2 + mom_z**2)
            particle_len = 200.

            end_vtx_r = self.cur_vtx_r*self.scale + mom_r/mom_abs*particle_len
            end_vtx_z = self.cur_vtx_z*self.scale + mom_z/mom_abs*particle_len

            dot_vtx = _get_dot_vertex(particle_barcode, end_vtx_r, end_vtx_z, is_dummy=True)
            self.dotfile.write(dot_vtx)

        particle_dot = _get_dot_particle(self.cur_vtx_barcode,
                                         end_vtx_barcode,
                                         particle_barcode,
                                         particle_id,
                                         particle_energy)
        self.dotfile.write(particle_dot)

    def close(self):
        self._end_opened_event()
        self.dotfile.close()

    def __del__(self):
        self.close()

    def _begin_event(self, raw_hepmc_line):
        self.event_open = True
        evt_num_column = 1
        evt_num = raw_hepmc_line.split()[evt_num_column]
        self.dotfile.write("digraph event_%s {\n" % evt_num)

    def _end_event(self):
        self.dotfile.write("}\n")

    def _end_opened_event(self):
        if self.event_open:
            self._end_event()
        self.event_open = False


def main(hepmc_file='events.hepmc', dot_file='events.dot'):
    begin_event_pattern = re.compile(r'^E .*$')
    vertex_pattern = re.compile(r'^V .*$')
    particle_pattern = re.compile(r'^P .*$')

    event_open = False
    with open(hepmc_file, 'r') as hepmc:
        dot = HepDotWriter(dot_file)

        for line in hepmc:
            if begin_event_pattern.match(line):
                dot.newEvent(line)
            elif vertex_pattern.match(line):
                dot.newVertex(line)
            elif particle_pattern.match(line):
                dot.addOutgoingParticle(line)
            else:
                pass # ignore unknown lines


if __name__ == '__main__':
    main()
