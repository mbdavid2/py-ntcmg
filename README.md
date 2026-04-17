# pyntcm

Python library wrapping the ntcm-g cpp implementation @ https://github.com/lguldur/ntcmg/

## Requirements

- g++
- python3-dev

## Setup

```bash
pip install .
```

## Changes

An additional `vtec()` method has been added along `stec()`. The method is the same but returns the vtec value directly instead of multiplying it by `vtec2stec_mapping_function(elev)`.

It takes only lat and lon (internally uses those for rx and sat and `satalt=20200000.0` ,`rxalt=0.0`).

## Example

```python
from pyntcm import _c_ext
result = _c_ext.stec(
            1.0, 2.0, 3.0,  # ai0, ai1, ai2
            0.5, 0.5, 0.0,  # rxlat, rxlon, rxalt
            0.6, 0.6, 20200000.0, # satlat, satlon, satalt
            12.5, 150.0     # utc_time, doy
        )
```

Both methods can also use keywords:

```python
from pyntcm import _c_ext
result = _c_ext.vtec(ai0=1.0, ai1=0.0, ai2=0.0,
                     lat=0.9, lon=0.2,
                     utc_time=12.0, doy=180.0)
```

# ntcmg

Implementation of the NTCM-G algorithm as described in the "EUROPEAN GNSS (GALILEO) OPEN SERVICE/NTCM-G IONOSPHERIC MODEL DESCRIPTION" document available here:

https://www.gsc-europa.eu/sites/default/files/NTCM-G_Ionospheric_Model_Description_-_v1.0.pdf

Forked from:

https://github.com/lguldur/ntcmg/
