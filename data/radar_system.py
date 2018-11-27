
class beam_pattern:
    # az and el are the pointing direction
    # az0 and el0 are the on-axis direction
    def __init__(self,gain_func,I_0,freq,beam_name=''):
        self.I_0=I_0
        self.freq=freq
        self.az0=None
        self.el0=None
        self.gain_func = gain_func
        self.beam_name = beam_name

    # move antenna towards direction az0 el0
    def point(self,az0,el0):
        self.az0=az0
        self.el0=el0

    # directly apply on-axis position
    def point_k(self,k):
        pass

    # get angle to on axis direction
    def angle(self,az,el):
        pass
    
    def gain_k(self,k):
        return self.gain_func(k,self)
        #return( airy(n.pi*coord.angle_deg(self.on_axis,k)/180.0,a=self.a0,f=self.f,I_0=self.I_0))

class antenna_array_beam(beam_pattern):
    def __init__(self,gain_func,antenna_positions,antenna_subgrouping,I_0,freq,beam_name=''):
        beam_pattern.__init__(self,gain_func,I_0,freq,beam_name)
        self.antenna_positions = antenna_positions
        self.antenna_subgrouping = antenna_subgrouping
        self.subgroup_positions = n.mean(antenna_positions[antenna_subgrouping]) #something like this

def subgroup_summation(k,beam):
    pass

def MU_radar_model():
    MU_beam = antenna_array_beam(subgroup_summation,
        antenna_positions, #this would be all antenna positions
        antenna_subgrouping, #this would be indexes for what antennas go into what subgroup
        I_0=I_0,
        freq=47e6,
        beam_name='MU-radar beam')

class rx_station:
    def __init__(self,name,lat,lon,alt,el_thresh,freq,rx_noise,beam):
        self.name = name
        self.lat = lat
        self.lon = lon
        self.el_thresh = el_thresh
        self.rx_noise = rx_noise
        self.beam = beam
        self.freq = freq
        self.wavelength = c.c/freq
        #self.ecef = coord.geodetic2ecef(lat,lon,0.0)
    
        self.__station_mode = 'RX'

    def __str__(self):
        string = self.__station_mode + ' '
        string+= "Antenna {}\n".format(self.name)
        string+= "{}:{}\n".format('Local horizon', str(el_thresh))
        string+= "{}:{}\n".format('Radiation pattern', str(ant))
        return(string)
    
class tx_station(rx_station):
    def __init__(self,name,lat,lon,alt,el_thresh,freq,rx_noise,ant,tx_power,tx_bandwidth,duty_cycle):
        rx_antenna.__init__(self,name,lat,lon,alt,el_thresh,freq,rx_noise,ant)
        self.__station_mode = 'TX'

        self.tx_bandwidth = tx_bandwidth
        self.duty_cycle = duty_cycle # max duty cycle
        self.tx_power = tx_power

class radar_system():
    def __init__(self,tx_lst,rx_lst,name):
      self.tx = tx_lst
      self.rx = rx_lst
      self.name = name

    def set_FOV(self,max_on_axis,horizon_elevation):
        self.max_on_axis=max_on_axis
        self.horizon_elevation = horizon_elevation
        for rx in self.rx:
            rx.el_thresh = horizon_elevation
        for tx in self._tx:
            tx.el_thresh = horizon_elevation

    def set_TX_bandwith(self,bw):
        for tx in self.tx:
            tx.tx_bandwidth = bw

    def set_beam(self,beam,mode=''):
        if mode=='TX':
            for tx in self.tx:
                tx.beam = beam
        elif mode=='RX':
            for rx in self.rx:
                rx.beam = beam
        else:
            for rx in self.rx:
                rx.beam = beam
            for tx in self.tx:
                tx.beam = beam

