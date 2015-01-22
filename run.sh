#!/bin/sh
exec R -e "shiny::runApp('.', launch.browser=F, port=8081)"
