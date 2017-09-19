#!/bin/bash

LL=$(stat -c %Y "$1")
while [ 1 == 1 ]; do
	[ -f "$1" ] || break
	LL2=$(stat -c %Y "$1")
	if [ $LL != $LL2 ]; then
		stat -c '%n %y' "$1"
		python "$1" || echo OOOPS
		LL=$LL2
	fi
	sleep 0.25
done
