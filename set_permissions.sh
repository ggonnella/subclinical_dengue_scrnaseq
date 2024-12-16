#!/bin/bash
find /srv/baia/inden -user $USER | xargs chmod g+rwX
find /srv/baia/inden -user $USER | xargs chgrp baia
