use p224::*;
use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() <= 1 || args.len() > 3 {
        println!("Usage: {} <private key> [public key]", args[0]);
        println!("Examples:");
        println!("  {} af0181e90420508e39e9d862f1680dc22f5f024fcfc11939d7c6fedb", args[0]);
        println!("  {} 14447aec7d78c691d0c1c94da1a6a85d9eefeddf8b42f51aa227376c", args[0]);
        println!("  {} af0181e90420508e39e9d862f1680dc22f5f024fcfc11939d7c6fedb 457e0910e4a34933c1ad3034ad92504c8324b8701e56c37716bf541967813d3ff1390e5c1d0c833f1fce3ff8bc69b277a072c3c31b239aacf", args[0]);
        println!("  {} 14447aec7d78c691d0c1c94da1a6a85d9eefeddf8b42f51aa227376c 4bcf74addac6c83a587eeb6d2a724158cebaecfed0af82a90434268e03c82f21e137c7341a70c0044ed058d5fe6c7aa38eb16542fdf5ae111", args[0]);
        return;
    }

    let key1 = Key::from_private_hex(&args[1]);
    println!("private key 1: {:?}", key1.get_private_hex());
    println!("public key 1:  {:?}", key1.get_public_hex());

    if args.len() == 3 {
        let key2 = Key::from_public_hex(&args[2]);
        let shared = key1.derive(&key2);
        println!("public key 2:  {:?}", key2.get_public_hex());
        println!("shared secret: {:?}", shared.get_public_hex());
    }
}
