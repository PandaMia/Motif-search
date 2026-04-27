# Deploy

Short deployment guide for `motif-search` on a DigitalOcean Droplet using Docker Compose, with Caddy running on the host.

## 1. Connect to the server

```bash
ssh deploy@YOUR_SERVER_IP
```

## 2. Go to the project directory

```bash
cd /srv/motif-search
```

## 3. Verify Docker is installed

Check:

```bash
docker --version
docker compose version
```

## 4. Start the application container

The application is published only to the host loopback interface:

```text
127.0.0.1:18080 -> container:8080
```

Start the service:

```bash
docker compose up -d --build
```

Check status:

```bash
docker compose ps
docker compose logs -f
```

Optional local check on the server:

```bash
curl -I http://127.0.0.1:18080
```

## 5. Configure host Caddy

Add this site block to the host Caddy configuration:

```caddy
motif.pandamia.org {
    encode zstd gzip
    reverse_proxy 127.0.0.1:18080
}
```

Then reload host Caddy:

```bash
sudo systemctl reload caddy
```

Caddy will automatically manage HTTPS as long as:
- `motif.pandamia.org` points to your server IP
- ports `80` and `443` are open

## 6. Update the service

```bash
cd /srv/motif-search
git pull
docker compose up -d --build
```

## 7. Useful commands

Stop the service:

```bash
docker compose down
```

Restart it:

```bash
docker compose restart
```

View logs:

```bash
docker compose logs -f
```

Check app reachability from the host:

```bash
curl -I http://127.0.0.1:18080
```
