SELECT COUNT(*) AS pickups
FROM dbo.nodes
WHERE instance_id = 4 AND node_type = 'Pickup Node';

SELECT MIN(node_id), MAX(node_id)
FROM dbo.nodes
WHERE instance_id = 4;