#!/usr/bin/env python

import math
import re # regex
import numpy as np


def _get_dot_vertex(vtx_id, r, z, is_dummy=False):
    """
    Generates a DOT formatted vertex for the given id and coordinates
    """
    if is_dummy:
        dummy = 'dummy_'
        attrib = 'shape=none,label=""'
    else:
        dummy = ''
        attrib = r'label="vtx {ID}\nr={r:.2f},z={z:.2f}"'.format(ID=vtx_id, r=r, z=z)

    dot = '    V_{dummy}{abs_id} [{attrib},pos="{zpos:.3f},{rpos:.3f}!"];\n'.format(
                                                                                dummy=dummy,
                                                                                abs_id=abs(int(vtx_id)),
                                                                                attrib=attrib,
                                                                                zpos=z,
                                                                                rpos=r)
    return dot


class HepDotWriter(object):
    """
    Generates a dot file representing the given particles, interaction vertices and events
    """

    def __init__(self, dotfile): #TODO: implement next:, vtx_threshold=np.nan, scale=1.):
        self.dotfile = open(dotfile, 'w')

        self.event_open = False
        self.cur_vtx = None
        self.cur_vtx_r = None
        self.cur_vtx_z = None

        # TODO: implement next:
        #self.vtx_threshold = vtx_threshold
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
        vtx_id_column = 1
        vtx_id = int(hepmc[vtx_id_column])
        vtx_abs_id = abs(vtx_id)

        # TODO: implement next
        #if vtx_abs_id > self.vtx_threshold:
        #    return

        x_column = 3
        x = float(hepmc[x_column])
        y_column = 4
        y = float(hepmc[y_column])
        z_column = 5

        self.cur_vtx_z = float(hepmc[z_column])
        self.cur_vtx_r = math.sqrt(x**2 + y**2)

        self.cur_vtx = 'V_{abs_id}'.format(abs_id=vtx_abs_id)
        dot_vtx = _get_dot_vertex(vtx_id,
                                  self.cur_vtx_r*self.scale,
                                  self.cur_vtx_z*self.scale)
        self.dotfile.write(dot_vtx)

    def addOutgoingParticle(self, raw_hepmc_line):
        line = raw_hepmc_line.split()

        particle_num_column = 1
        particle_num = line[particle_num_column]

        mom_x_column = 3
        mom_x = float(line[mom_x_column])
        mom_y_column = 4
        mom_y = float(line[mom_y_column])
        mom_z_column = 5
        mom_z = float(line[mom_z_column])

        end_vtx_num_column = 11
        end_vtx_num = abs(int(line[end_vtx_num_column]))

        # TODO: implement next
        #if end_vtx_num > self.vtx_threshold:
        #    return

        if not end_vtx_num:
            # create dummy end node for partiles that don't have end vertices

            mom_r = math.sqrt(mom_x**2 + mom_y**2)
            mom_abs = math.sqrt(mom_r**2 + mom_z**2)
            particle_len = 200.

            end_vtx_r = self.cur_vtx_r*self.scale + mom_r/mom_abs*particle_len
            end_vtx_z = self.cur_vtx_z*self.scale + mom_z/mom_abs*particle_len

            dot_vtx = _get_dot_vertex(particle_num, end_vtx_r, end_vtx_z, is_dummy=True)
            self.dotfile.write(dot_vtx)

            end_vtx = 'V_dummy_{num}'.format(num=particle_num)
        else:
            end_vtx = 'V_{num}'.format(num=end_vtx_num)

        particle_dot = '    {cur_vtx} -> {end_vtx} [label="p {num}"];\n'.format(
                                                                            cur_vtx=self.cur_vtx,
                                                                            end_vtx=end_vtx,
                                                                            num=particle_num)
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
