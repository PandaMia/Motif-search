# Deploy

Short deployment guide for `motif-search` on a DigitalOcean Droplet using Docker Compose and Caddy.

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

## 4. First launch over HTTP

The current `Caddyfile` is already configured for plain HTTP on `:80`, so no changes are needed for the first launch.

Start the service:

```bash
docker compose up -d --build
```

Check status:

```bash
docker compose ps
docker compose logs -f
```

After that, the service should be available by server IP:

```text
http://YOUR_SERVER_IP
```

## 5. Switch to a domain and HTTPS

1. Create an `A` record and point it to your server IP.
2. Open `Caddyfile`.
3. Replace this line:

```caddy
:80 {
```

with:

```caddy
your-domain.com {
```

or:

```caddy
motif.example.com {
```

4. Reload the containers:

```bash
docker compose up -d
```

Caddy will automatically issue and renew the TLS certificate.

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
