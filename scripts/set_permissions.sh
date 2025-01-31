#!/bin/bash
find /srv/baia/prj/subclinical_dengue_scrnaseq -user $USER | xargs chmod g+rwX
find /srv/baia/prj/subclinical_dengue_scrnaseq -user $USER | xargs chgrp baia
