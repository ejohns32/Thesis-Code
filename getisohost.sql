SELECT
    iso.isoID,
    host.commonName
    FROM
        Isolates as iso,
        Host as host 
    WHERE
        iso.hostID = host.hostID
;
